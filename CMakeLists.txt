cmake_minimum_required(VERSION 3.9)

project(TWSM LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache) # Less useful to do it for linking, see edit2
endif(CCACHE_FOUND)

find_package(NetCDF REQUIRED)
find_package(OpenMP)
# use -fsanitize=address (above) to check for memory errors


add_subdirectory(common/richdem EXCLUDE_FROM_ALL)

add_executable(twsm.x
  Combined/ArrayPack.cpp
  Combined/parameters.cpp
  Combined/transient_groundwater.cpp
  Combined/TWSM.cpp
)
target_link_libraries(twsm.x PRIVATE richdem ${NETCDF_LIBRARIES} OpenMP::OpenMP_CXX)
target_compile_features(twsm.x PRIVATE cxx_std_17)
target_compile_options(twsm.x PRIVATE -Wall -pedantic)

add_executable(dephier.x
  Combined/run_dephier.cpp
  Combined/parameters.cpp
  Combined/ArrayPack.cpp
)
target_link_libraries(dephier.x PRIVATE richdem ${NETCDF_LIBRARIES} OpenMP::OpenMP_CXX)
target_compile_features(dephier.x PRIVATE cxx_std_17)
target_compile_options(dephier.x PRIVATE -Wall -pedantic)
