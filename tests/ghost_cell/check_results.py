#!/usr/bin/env python3
"""
Compare 1-process and 2-process WTM output TIFs.

With the ghost-cell fix both runs must agree at all interior cells.
Without it the two halves of the domain are hydrologically decoupled at the
MPI processor boundary (x = NX//2), producing a clearly visible (~O(1) m) error.
"""

import sys
import os
import glob
import numpy as np
import rasterio

NX = 12


def last_tif(prefix, outdir):
    pattern = os.path.join(outdir, f"{prefix}*.tif")
    tifs = sorted(glob.glob(pattern))
    if not tifs:
        raise FileNotFoundError(f"No TIF files matching {pattern}")
    return tifs[-1]


def load(path):
    with rasterio.open(path) as src:
        data = src.read(1).astype(np.float64)
    return data


def main():
    dir1p = "out_1p"
    dir2p = "out_2p"

    f1 = last_tif("out_", dir1p)
    f2 = last_tif("out_", dir2p)
    print(f"1-process output : {f1}")
    print(f"2-process output : {f2}")

    w1 = load(f1)
    w2 = load(f2)

    if w1.shape != w2.shape:
        print(f"FAIL: shape mismatch {w1.shape} vs {w2.shape}", file=sys.stderr)
        sys.exit(1)

    diff = np.abs(w1 - w2)

    # Interior land cells only (exclude ocean edges at row 0, row NY-1, col 0, col NX-1)
    interior = diff[1:-1, 1:-1]

    max_diff      = interior.max()
    mean_diff     = interior.mean()
    boundary_diff = diff[1:-1, NX // 2 - 1 : NX // 2 + 1].max()

    print(f"Max  |Δwtd| interior     : {max_diff:.6f} m")
    print(f"Mean |Δwtd| interior     : {mean_diff:.6f} m")
    print(f"Max  |Δwtd| at MPI bound : {boundary_diff:.6f} m")

    # Threshold: numerical summation-order differences are O(1e-10) m.
    # The ghost-cell bug produces O(1) m errors near the boundary.
    TOLERANCE = 1e-4  # generous; anything above ~0.01 m signals the bug

    if max_diff > TOLERANCE:
        print(
            f"\nFAIL: max difference {max_diff:.4f} m exceeds tolerance {TOLERANCE} m.\n"
            "Ghost-cell error is present — the MPI boundary suppresses inter-rank flux.",
            file=sys.stderr,
        )
        sys.exit(1)
    else:
        print(f"\nPASS: 1-process and 2-process outputs agree to within {TOLERANCE} m.")
        sys.exit(0)


if __name__ == "__main__":
    main()
