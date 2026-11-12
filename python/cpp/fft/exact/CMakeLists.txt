# The exact FFT add-on: FFTW 3.3.8 (GPL-2.0-or-later) with MATLAB R2024b's codelet sets, as two shared libraries
# (single and double precision) that the port loads at run time for MATLAB's bits (cpp/fft/backend.cpp).
#
#   as its own package:   pip install python/cpp/fft/exact      (petakit5d-fftw; libraries in petakit5d_fftw/)
#   in the port's build:  cmake -S python -B <build> -DPK5D_FFT_EXACT=ON   (libraries in <build>/fftw_exact;
#                         PETAKIT5D_FFTW=<build>/fftw_exact)
# Neither installs anything into the petakit5d package: the port as distributed contains no GPL code.
cmake_minimum_required(VERSION 3.26...4.2)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(petakit5d_fftw LANGUAGES C)
    if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
        set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    endif()
    set(_pk5d_fftw_standalone TRUE)
    set(PK5D_FFTW_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/lib)
else()
    enable_language(C)
    set(_pk5d_fftw_standalone FALSE)
    set(PK5D_FFTW_OUTPUT_DIR ${CMAKE_BINARY_DIR}/fftw_exact)
endif()

include(${CMAKE_CURRENT_SOURCE_DIR}/fftw.cmake)

if(_pk5d_fftw_standalone AND SKBUILD)
    # the add-on's wheel: the libraries next to petakit5d_fftw/__init__.py, FFTW's license with them
    install(TARGETS pk5d_fftw3f pk5d_fftw3 LIBRARY DESTINATION petakit5d_fftw RUNTIME DESTINATION petakit5d_fftw)
    install(FILES ${pk5d_fftw_SOURCE_DIR}/COPYING DESTINATION petakit5d_fftw RENAME FFTW-COPYING)
endif()
