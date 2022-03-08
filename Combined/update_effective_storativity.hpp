#include "ArrayPack.hpp"
#include "parameters.hpp"

//Deal with the fact that porosity is 1 above ground but [value] below ground. It is included directly in the matrix calculation,
//so we can't just scale water change before or after. Instead, we are scaling the effective storativity according to how much
//change in water table we are projecting to see during that time step.
double updateEffectiveStorativity(const double my_original_wtd, const double my_wtd_T, const double my_porosity, const double starting_effective_storativity){
  double projected_full_step_wtd = my_original_wtd + (my_wtd_T - my_original_wtd)*2.;
  if(my_original_wtd <= 0. && projected_full_step_wtd <= 0.) //both are below ground, so we can use the original porosity
    return static_cast<double>(my_porosity);
  else if(my_original_wtd >= 0. && projected_full_step_wtd >= 0.) //both are above ground, so the porosity is 1
    return 1.;
  else{
    double change_in_water_column_thickness = std::abs(projected_full_step_wtd - my_original_wtd);
    if(my_original_wtd < 0. && projected_full_step_wtd > 0.){ //started below ground and ended above ground
      // First, scale the change in wtd as if the whole column has the porosity of the belowground area
      double scaled_change_in_water_column_thickness = change_in_water_column_thickness * starting_effective_storativity/static_cast<double>(my_porosity);
      // Then get just that portion that is >0 (aboveground), and scale it back down to the equivalent thickness with 100% porosity (surface water)
      double aboveground_water_column_thickness = (scaled_change_in_water_column_thickness + my_original_wtd) * starting_effective_storativity;//Above-ground effective porosity is ==1, so no need to actually divide by 1 here (save on computation).
      //belowground water column thickness is equal to -original_wtd, so no need to assign it to a new variable.
      return (aboveground_water_column_thickness + static_cast<double>(my_porosity)*-my_original_wtd) / (aboveground_water_column_thickness - my_original_wtd);
    }
    else{ //This means that (my_original_wtd > 0. && projected_full_step_wtd < 0.). started above ground and ended below ground.
      double scaled_change_in_water_column_thickness = change_in_water_column_thickness * starting_effective_storativity; //divided by 1 for above-ground porosity
      if(scaled_change_in_water_column_thickness < projected_full_step_wtd){  //when rescaling the water according to the porosity values, there is no longer enough to reach below ground; so the scaled new porosity will be = 1
        return 1;
      }
      else{
        // Get the belowground water thickness and expand it to a deeper depth in order to account for changing porosity
        double belowground_water_column_thickness = (scaled_change_in_water_column_thickness - my_original_wtd) * starting_effective_storativity/static_cast<double>(my_porosity);
        return (my_original_wtd + static_cast<double>(my_porosity)*belowground_water_column_thickness) / (my_original_wtd + belowground_water_column_thickness);
      }
    }
  }
}