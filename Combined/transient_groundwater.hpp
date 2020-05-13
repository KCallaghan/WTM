#ifndef transient_groundwater_h
#define transient_groundwater_h

#include <stdint.h>
#include "../common/netcdf.hpp"
#include "ArrayPack.hpp"
#include "parameters.hpp"
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <richdem/common/Array2D.hpp>
#include <richdem/common/timer.hpp>
#include <richdem/common/ProgressBar.hpp>
#include <stdexcept>
#include <string>
#include <vector>
#include <fstream>

using namespace std;

typedef std::vector<double> dvec;
typedef rd::Array2D<float>  f2d;

class FanDarcyGroundwater{

public:
    ////////////////////////////////////
    // INSTANCE VARIABLES AND OBJECTS //
    ////////////////////////////////////

    //char infname[] = "path";
    //char *argv[2] = {&infname, &infname}

    ArrayPack arp;
    Parameters params;

    /////////////////
    // CONSTRUCTOR //
    /////////////////

    FanDarcyGroundwater();
    FanDarcyGroundwater(Parameters _params, ArrayPack _arp);

    ///////////////
    // FUNCTIONS //
    ///////////////
    void set_arp(ArrayPack &_arp);
    void set_params(Parameters &_params);

    void initialize();
    void update(bool _log=true);
    void run();
    void finalize();


private:
    //////////////////////
    // INSTANCE OBJECTS //
    //////////////////////



    ////////////////////////
    // INSTANCE VARIABLES //
    ////////////////////////

    double transmissivityN;
    double transmissivityS;
    double transmissivityW;
    double transmissivityE;

    double headCenter;
    double headN;
    double headS;
    double headW;
    double headE;

    double wtdCenter;
    double wtdN;
    double wtdS;
    double wtdW;
    double wtdE;

    // These variables are used to monitor the state of the calculation
    double total_changes  = 0.;
    float max_total       = 0.;
    float min_total       = 0.;
    float max_change      = 0.;


    ///////////////
    // FUNCTIONS //
    ///////////////

    /**
     * @brief Provides transmissivity within each cell [m^2/s].
     * @details Transmissivity [units = m^2/s], which is the hydraulic
     *          conductivity, integrated vertically from deep depth up to the
     *          water table. It changes through time as the water-table
     *          depth changes. This is due to a combination of:
     *          (a) the increasing to total thickness of the water column, and
     *          (b) the exponentially-increasing hydraulic conductivity that
     *              reaches its maximum value at 1.5 meters depth, at which
     *              point it becomes a constant, based on soil-survey data.
     *          We follow Ying Fan Reinfelder et al., 2013, Science,
     *          "Global patterns of groundwater table depth", in our approach.
     *          This combines the aforementioned expoenential decay in
     *          hydraulic conductivity below a mapped "soil" layer.
     * @param x The x-coordinate of the cell in question
     * @param y The y-coordinate of the cell in question
     * @param ArrayPack Global arrays. Here we use:
     *       - fdepth: The e-folding depth based on slope and temperature.
     *                 This describes the decay of kcell with depth.
     *       - wtd:    The water-table depth. We use a different calculation
     *                 for water tables above vs below 1.5 m depth.
     *       - ksat:   Surface hydraulic conductivity, based on soil types
     * @return  The transmissivity value for the cell in question. This is the
     *          integration of the hydraulic conductivity from -infinity to the
     *          groundwater table.
     */
    double computeTransmissivity(uint32_t x, uint32_t y);

    /**
     * @brief Compute mean transmissivity between neighboring cells (N, S, W, E)
     */
    void computeNeighborTransmissivity(uint32_t x, uint32_t y);

    /**
     * @brief Returns the maximum value in an array (max size 256 items)
     */
    double computeArrayMax(double *_val[], uint8_t size);

    /**
     * @brief Returns the minimum value in an array (max size 256 items)
     */
    float computeArrayMin(float *_val[], uint8_t size);

    /**
     * @brief Returns the maximum stable time step with a 2x factor of safety
     * @details Uses a 2D diffusion von Neumann stability analysis using the
     *          "worst case" scenario highest transmissivity, combined with
     *          a porosity-based amplification factor.
     */
    double computeMaxStableTimeStep(uint32_t x, uint32_t y);

    /**
     * @brief Calculates water-table depth change at a cell (and associated
     *        surrounding cells) and updates the class variables associated
     *        with these.
     */
    void computeWTDchangeAtCell(int32_t x, int32_t y, double dt);

    /**
     * @brief Updates the wtd_depth_total array at a cell(x,y) using the
     * pre-set time step and dynamic time stepping within this as needed.
     */
    void updateCell(uint32_t x, uint32_t y);

    /**
     * @brief Create a continuously updating output file with basic information
     *        on the run.
     */
    void logToFile();
};

#endif
