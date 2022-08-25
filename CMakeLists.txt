cmake_minimum_required(VERSION 3.9)

project(WTM LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# use -fsanitize=address to check for memory errors
set(CMAKE_CXX_FLAGS_RELTESTING "-O3 -g -Wall -Wextra -pedantic -Wshadow")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -Wall -Wextra -pedantic -Wshadow")

# RichDEM requires GDAL
set(USE_GDAL ON)

# # TODO: Enable to perform many linting checks
# set(CMAKE_CXX_CLANG_TIDY
#         clang-tidy;
#             # -header-filter=.*;
#             -checks=*,-llvmlibc-restrict-system-libc-headers,-modernize-use-trailing-return-type,-llvmlibc-callee-namespace,-bugprone-branch-clone,-google-readability-namespace-comments,-llvm-namespace-comment,-llvmlibc-implementation-in-namespace,-readability-uppercase-literal-suffix,-readability-function-cognitive-complexity,-fuchsia-default-arguments-calls,-cppcoreguidelines-pro-bounds-array-to-pointer-decay,-hicpp-uppercase-literal-suffix,-hicpp-no-array-decay,-misc-unused-parameters;
#             -warnings-as-errors=*;)

find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache) # Less useful to do it for linking, see edit2
endif(CCACHE_FOUND)

find_package(GDAL REQUIRED)
#find_package(MPI REQUIRED)
find_package(OpenMP)
find_package(PkgConfig REQUIRED)
find_package(Sanitizers)

# NOTE: You may need to update this path to identify PETSc's location
set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:/usr/lib/petscdir/petsc3.17.1/lib/pkgconfig/")
pkg_check_modules(PETSC PETSc>=3.17.1 IMPORTED_TARGET REQUIRED)
message(STATUS "Found PETSc ${PETSC_VERSION}")

add_subdirectory(common/richdem EXCLUDE_FROM_ALL)
add_subdirectory(common/fmt EXCLUDE_FROM_ALL)

# Define a library aggregating all the model componenets that we can use for
# both running simulations and running tests
add_library(wtm
  src/ArrayPack.cpp
  src/irf.cpp
  src/parameters.cpp
  src/transient_groundwater.cpp
  src/update_effective_storativity.cpp
)
target_compile_options(wtm PRIVATE -Wfloat-conversion -Wall -Wextra -pedantic -Wshadow)
target_link_libraries(wtm PUBLIC richdem OpenMP::OpenMP_CXX fmt::fmt PkgConfig::PETSC)
target_compile_features(wtm PUBLIC cxx_std_20)

add_executable(wtm.x
  src/WTM.cpp
)
target_compile_options(wtm.x PRIVATE -Wfloat-conversion -Wall -Wextra -pedantic -Wshadow)
target_link_libraries(wtm.x PRIVATE wtm)
target_compile_features(wtm.x PRIVATE cxx_std_20)
target_compile_definitions(wtm PRIVATE DOCTEST_CONFIG_DISABLE)

add_sanitizers(wtm.x)

add_executable(dephier.x
  src/run_dephier.cpp
  src/parameters.cpp
  src/ArrayPack.cpp
)
target_compile_options(dephier.x PRIVATE -Wfloat-conversion -Wall -Wextra -pedantic -Wshadow)
target_link_libraries(dephier.x PRIVATE richdem OpenMP::OpenMP_CXX fmt::fmt)
target_compile_features(dephier.x PRIVATE cxx_std_20)
add_sanitizers(dephier.x)

add_executable(tests.x
  src/parameters.cpp
  src/tests.cpp
)
target_compile_options(tests.x PRIVATE -Wfloat-conversion -Wall -Wextra -pedantic -Wshadow)
target_link_libraries(tests.x PRIVATE wtm)
target_compile_features(tests.x PRIVATE cxx_std_20)
