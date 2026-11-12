# Build of the petakit5d extension modules: one nanobind module per kernel group,
# petakit5d._ext.<group>, from python/cpp/<group>/ (see cmake/Pk5dExtension.cmake and
# docs/decisions/0001-build.md).
#
# Two ways to build:
#   pip install .                         scikit-build-core; modules go into the wheel
#   cmake -S python -B <build> -G Ninja   then cmake --build <build>; modules land in
#                                         <build>/petakit5d/_ext, used with
#                                         PYTHONPATH=python/src PETAKIT5D_EXT_DIR=<build>/petakit5d/_ext
cmake_minimum_required(VERSION 3.26...4.2)

project(petakit5d LANGUAGES C CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(PK5D_REQUIRE_OPENMP "Fail when no OpenMP is found (OFF: build serial kernels)" ON)
set(PK5D_EXT_OUTPUT_DIR "${CMAKE_BINARY_DIR}/petakit5d/_ext" CACHE PATH
    "Where the build tree puts the extension modules (PETAKIT5D_EXT_DIR)")

find_package(Python 3.11 REQUIRED COMPONENTS Interpreter Development.Module
             OPTIONAL_COMPONENTS Development.SABIModule)

if(NOT nanobind_ROOT)
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT
    RESULT_VARIABLE _pk5d_nb_res)
  if(NOT _pk5d_nb_res EQUAL 0)
    message(FATAL_ERROR "nanobind not found for ${Python_EXECUTABLE} (pip install nanobind)")
  endif()
endif()
find_package(nanobind CONFIG REQUIRED)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(Pk5dOpenMP)
include(Pk5dExtension)

add_subdirectory(cpp)

pk5d_finalize()
