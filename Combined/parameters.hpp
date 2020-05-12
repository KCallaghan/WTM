#ifndef _parameters_hpp_
#define _parameters_hpp_

#include <cmath>
#include <limits>
#include <string>

const std::string UNINIT_STR = "uninitialized";

class Parameters{
public:
    Parameters(const std::string config_file);

    int32_t iterations = -1;
    int32_t maxiter    = -1;

    std::string surfdatadir  = UNINIT_STR;
    std::string region       = UNINIT_STR;
    std::string run_type     = UNINIT_STR;
    std::string time_start   = UNINIT_STR;
    std::string time_end     = UNINIT_STR;
    std::string textfilename = UNINIT_STR;
    std::string outfilename  = UNINIT_STR;

    int16_t cells_per_degree = -1;

    const float64_t UNDEF  = -1.0e7;

    bool infiltration_on;

    float64_t southern_edge         = std::numeric_limits<float64_t>\
                                         ::signaling_NaN();
    float64_t deltat                = std::numeric_limits<float64_t>\
                                         ::signaling_NaN();
    float64_t cellsize_n_s_metres   = std::numeric_limits<float64_t>\
                                         ::signaling_NaN();
    float32_t  infiltration         = 0.;
    int32_t    cycles_done          = 0;
    float32_t  total_wtd_change     = 0.;
    float32_t  wtd_mid_change       = 0.;
    float32_t  GW_wtd_change        = 0.;
    float32_t  abs_total_wtd_change = 0.;
    float32_t  abs_wtd_mid_change   = 0.;
    float32_t  abs_GW_wtd_change    = 0.;
    float32_t  infiltration_change  = 0.;
    float32_t  surface_change       = 0.;
    int32_t    total_cycles         = -1;

    //Set for convenience within the code
    int32_t ncells_x  = -1;
    int32_t ncells_y  = -1;

    void print() const;
};

#endif
