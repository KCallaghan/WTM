cmake_minimum_required(VERSION 3.9)

project(WTM LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(CMAKE_CXX_FLAGS_RELTESTING "-O3 -g -Wall -Wextra -pedantic -Wshadow -lpetsc")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -Wall -Wextra -pedantic -Wshadow -lpetsc")

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

find_package(PkgConfig REQUIRED)
find_package(OpenMP)
find_package(Sanitizers)
find_package(GDAL REQUIRED)

pkg_search_module(PETSC PETSc>3.5.0 IMPORTED_TARGET REQUIRED)
# use -fsanitize=address (above) to check for memory errors

add_subdirectory(common/richdem EXCLUDE_FROM_ALL)
add_subdirectory(common/fmt EXCLUDE_FROM_ALL)

add_executable(wtm.x
  src/ArrayPack.cpp
  src/irf.cpp
  src/parameters.cpp
  src/transient_groundwater.cpp
  src/update_effective_storativity.cpp
  src/WTM.cpp
)
target_compile_options(wtm.x PRIVATE -Wfloat-conversion)
target_link_libraries(wtm.x PRIVATE ${PETSC_LIBRARIES} ${PETSC_PACKAGE_LIBS} richdem OpenMP::OpenMP_CXX fmt::fmt)
target_compile_features(wtm.x PRIVATE cxx_std_17)
target_compile_definitions(wtm.x PRIVATE DOCTEST_CONFIG_DISABLE)
add_sanitizers(wtm.x)

add_executable(dephier.x
  src/run_dephier.cpp
  src/parameters.cpp
  src/ArrayPack.cpp
)
target_link_libraries(dephier.x PRIVATE richdem OpenMP::OpenMP_CXX fmt::fmt)
target_compile_features(dephier.x PRIVATE cxx_std_17)
add_sanitizers(dephier.x)

add_executable(tests.x
  src/parameters.cpp
  src/test_fill_spill_merge.cpp
  src/tests.cpp
  src/transient_groundwater.cpp
  src/update_effective_storativity.cpp
)
target_link_libraries(tests.x PRIVATE richdem OpenMP::OpenMP_CXX fmt::fmt)
target_compile_features(tests.x PRIVATE cxx_std_17)
