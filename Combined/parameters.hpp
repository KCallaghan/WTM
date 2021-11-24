#ifndef _parameters_hpp_
#define _parameters_hpp_

#include <stdint.h>
#include <cmath>
#include <limits>
#include <string>

const std::string UNINIT_STR = "uninitialized";

class Parameters{
public:
    Parameters();
    Parameters(const std::string config_file);

    int32_t iterations = -1;
    int32_t maxiter    = -1;

    std::string outfilename  = UNINIT_STR;
    std::string region       = UNINIT_STR;
    std::string run_type     = UNINIT_STR;
    std::string surfdatadir  = UNINIT_STR;
    std::string textfilename = UNINIT_STR;
    std::string time_start   = UNINIT_STR;
    std::string time_end     = UNINIT_STR;

    int16_t cells_per_degree = -1;

    double UNDEF  = -1.0e7;

    bool infiltration_on;
    bool supplied_wt;
    bool evap_mode            = 1; // Default potential evaporation
    bool fsm_on               = 1; // Default surface water on

    double    deltat               = std::numeric_limits<double>::signaling_NaN();
    double    fdepth_a             = 0.;
    double    fdepth_b             = 0.;
    double    fdepth_fmin          = 0.;
    double    southern_edge        = std::numeric_limits<double>::signaling_NaN();
    int32_t   picard_iterations    = 1; //Default only one iteration
    int32_t   total_cycles         = -1;

    double    abs_GW_wtd_change    = 0.;
    double    abs_total_wtd_change = 0.;
    double    abs_wtd_mid_change   = 0.;
    double    cellsize_n_s_metres  = std::numeric_limits<double>::signaling_NaN();
    int32_t   cycles_done          = 0;
    double    total_wtd_change     = 0.;
    double    wtd_mid_change       = 0.;
    double    GW_wtd_change        = 0.;
    double    infiltration_change  = 0.;




      float s_big = 500;
      float s_medium = 10;
      float s_small = 2;
      float s_tiny = 1;
      float s_itsy = 0.2;
      float s_bitsy = 0.05;
      float s_between = 0.01;
      float s_spider = 0.00005;


    //Set for convenience within the code
    int32_t ncells_x  = -1;
    int32_t ncells_y  = -1;

    void print() const;
};

#endif
