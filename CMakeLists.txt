cmake_minimum_required(VERSION 3.9)

project(TWSM LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(CMAKE_CXX_FLAGS_RELTESTING "-O3 -g -Wall -pedantic")

find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache) # Less useful to do it for linking, see edit2
endif(CCACHE_FOUND)

find_package(OpenMP)
find_package(Sanitizers)
find_package(GDAL)
# use -fsanitize=address (above) to check for memory errors


add_subdirectory(common/richdem EXCLUDE_FROM_ALL)




if(USE_GDAL AND GDAL_FOUND)
  add_executable(twsm.x
    Combined/ArrayPack.cpp
    Combined/parameters.cpp
    Combined/transient_groundwater.cpp
    Combined/TWSM.cpp
  )
  target_link_libraries(twsm.x PRIVATE richdem OpenMP::OpenMP_CXX)
  target_compile_features(twsm.x PRIVATE cxx_std_17)
  target_compile_definitions(twsm.x PRIVATE DOCTEST_CONFIG_DISABLE)
  add_sanitizers(twsm.x)


  add_executable(dephier.x
    Combined/run_dephier.cpp
    Combined/parameters.cpp
    Combined/ArrayPack.cpp
  )
  target_link_libraries(dephier.x PRIVATE richdem OpenMP::OpenMP_CXX)
  target_compile_features(dephier.x PRIVATE cxx_std_17)
  add_sanitizers(dephier.x)
else()
  message(WARNING "GDAL is needed to build some of the FSM performance tests and example programs; it is not needed to use the FSM library in another project.")
  message(WARNING "Without GDAL, the WTM *will not build*.")
endif()

