#!/usr/bin/env python3
"""
Generate synthetic GeoTIFF inputs for the ghost-cell MPI validation test.

Domain: 102 x 3 cells (100 interior x-cells, 1 interior y-row).
Left half (x=1..50):  ksat = 1e-4 m/s
Right half (x=51..100): ksat = 1e-3 m/s  (10x higher)

The ksat discontinuity at x=50/51 coincides with the MPI processor boundary
when the job is split with -da_processors_x 2 -da_processors_y 1.  With the
ghost-cell bug the two halves are hydrologically decoupled (no flux across the
boundary); with the fix, the 2-process result matches the 1-process result.
"""

import numpy as np
import os
import rasterio
from rasterio.transform import from_bounds

NX = 12    # total columns (10 interior + 2 ocean-edge columns)
NY = 3     # total rows (1 interior row + 2 ocean-edge rows)

REGION = "ghost_cell_test"
TIME   = "t0"
OUTDIR = os.path.join(os.path.dirname(__file__), "inputs")
os.makedirs(OUTDIR, exist_ok=True)

# Arbitrary geographic extent; richdem only needs dimensions + data values.
transform = from_bounds(0, 0, NX, NY, NX, NY)
CRS = "EPSG:4326"


def write_tif(path, data, dtype="float32"):
    with rasterio.open(
        path, "w",
        driver="GTiff",
        height=NY, width=NX,
        count=1, dtype=dtype,
        crs=CRS,
        transform=transform,
    ) as dst:
        dst.write(data.astype(dtype), 1)


# Flat topography at 100 m.  Edge cells will be forced to 0 (ocean) by the
# model's land_mask.setEdges(0) + wtd=topo=0 for mask==0 logic.
topo = np.full((NY, NX), 100.0, dtype=np.float32)

# Flat slope → fdepth = fdepth_a (set to 200 in the config)
slope = np.zeros((NY, NX), dtype=np.float32)

# Land mask: 1 inside, 0 at edges.  The model also calls setEdges(0) itself,
# so this is just for the file to have valid dimensions.
mask = np.ones((NY, NX), dtype=np.float32)
mask[0, :] = 0
mask[-1, :] = 0
mask[:, 0] = 0
mask[:, -1] = 0

precip          = np.full((NY, NX), 0.3,  dtype=np.float32)   # m/yr
evap            = np.zeros((NY, NX),      dtype=np.float32)   # m/yr
open_water_evap = np.full((NY, NX), 0.4,  dtype=np.float32)   # m/yr
winter_temp     = np.zeros((NY, NX),      dtype=np.float32)   # deg C  (>-5 → no permafrost)

# Heterogeneous ksat: factor-of-10 jump at the midpoint.
# This asymmetry means there is non-zero flux at the MPI boundary in the
# correct solution; the ghost-cell bug suppresses it, causing a clear error.
ksat = np.full((NY, NX), 1e-4, dtype=np.float32)
ksat[:, NX // 2:] = 1e-3   # right half: 10× higher

porosity = np.full((NY, NX), 0.25, dtype=np.float32)

# Write all files
files = {
    f"{REGION}_{TIME}_topography.tif":            topo,
    f"{REGION}_{TIME}_slope.tif":                 slope,
    f"{REGION}_{TIME}_mask.tif":                  mask,
    f"{REGION}_{TIME}_precipitation.tif":         precip,
    f"{REGION}_{TIME}_evaporation.tif":           evap,
    f"{REGION}_{TIME}_open_water_evaporation.tif": open_water_evap,
    f"{REGION}_{TIME}_winter_temperature.tif":    winter_temp,
    f"{REGION}_horizontal_ksat.tif":              ksat,
    f"{REGION}_porosity.tif":                     porosity,
}

for fname, arr in files.items():
    path = os.path.join(OUTDIR, fname)
    write_tif(path, arr)
    print(f"  wrote {path}")

print("Done.")
