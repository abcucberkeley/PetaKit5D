# The exact FFT add-on: FFTW 3.3.8 (GPL-2.0-or-later) with MATLAB R2024b's codelet sets, as two shared libraries
# (single and double precision) that the port loads at run time for MATLAB's bits (cpp/fft/backend.cpp).
#
#   as its own package:   pip install python/cpp/fft/exact      (petakit5d-fftw; libraries in petakit5d_fftw/)
#   in the port's build:  cmake -S python -B <build> -DPK5D_FFT_EXACT=ON   (libraries in <build>/fftw_exact;
#                         PETAKIT5D_FFTW=<build>/fftw_exact)
# Neither installs anything into the petakit5d package: the port as distributed contains no GPL code.
cmake_minimum_required(VERSION 3.26...4.2)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    # Windows: the compilers of a source build, as petakit5d's, chosen before project() (cmake/Pk5dToolchain.cmake,
    # petakit5d's file: MinGW-builds' GCC 12.2.0 with the UCRT where it is found, else MSVC)
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Pk5dToolchain.cmake)
    if(CMAKE_HOST_WIN32)
        pk5d_windows_toolchain()
    endif()
    project(petakit5d_fftw LANGUAGES C)
    if(WIN32)
        message(STATUS "petakit5d-fftw: the Windows toolchain: ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION} "
                       "(${CMAKE_C_COMPILER}): ${PK5D_TOOLCHAIN_WHY}")
    endif()
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
    # MinGW-w64: the toolchain's notices of what the libraries link in (GCC's libgcc, MinGW-w64's runtime), copied
    # as the toolchain ships them into petakit5d_fftw/licenses/; with MinGW-builds' GCC 12.2.0 with the UCRT
    # (cmake/Pk5dToolchain.cmake's check; petakit5d's tools/get_mingw.sh) their SHA-256 must be that toolchain's
    # (petakit5d's cmake/Pk5dExtension.cmake does the same)
    if(WIN32 AND CMAKE_C_COMPILER_ID STREQUAL "GNU")
        get_filename_component(_pk5d_tc "${CMAKE_C_COMPILER}" DIRECTORY)
        get_filename_component(_pk5d_tc "${_pk5d_tc}/../licenses" ABSOLUTE)
        set(PK5D_MINGW_NOTICES_DIR "${_pk5d_tc}" CACHE PATH "MinGW-w64: the toolchain's licenses folder")
        pk5d_mingw_toolchain_check("${CMAKE_C_COMPILER}" _pk5d_pinned _pk5d_toolchain)
        foreach(_e
                "gcc/COPYING3=8ceb4b9ee5adedde47b31e975c1d90c73ad27b6b165a1dcd80c7c545eb65b903"
                "gcc/COPYING.RUNTIME=9d6b43ce4d8de0c878bf16b54d8e7a10d9bd42b75178153e3af6a815bdc90f74"
                "mingw-w64/COPYING.MinGW-w64-runtime.txt=36a4857b6557a97baa82762444ffd6a74c7575015595dc26ff02215c299f81f6")
            string(REPLACE "=" ";" _kv "${_e}")
            list(GET _kv 0 _rel)
            list(GET _kv 1 _sha)
            set(_f "${PK5D_MINGW_NOTICES_DIR}/${_rel}")
            if(NOT EXISTS "${_f}")
                message(WARNING "MinGW-w64: no ${_f}: the wheel would lack that notice (PK5D_MINGW_NOTICES_DIR)")
                continue()
            endif()
            file(SHA256 "${_f}" _got)
            if(_pk5d_pinned AND NOT _got STREQUAL _sha)
                message(FATAL_ERROR "MinGW-w64: ${_f} has SHA-256 ${_got}, not MinGW-builds GCC 12.2.0's ${_sha}")
            elseif(NOT _pk5d_pinned)
                message(STATUS "petakit5d-fftw MinGW-w64 notice ${_rel}: SHA-256 ${_got} (not MinGW-builds 12.2.0's)")
            endif()
            get_filename_component(_sub "${_rel}" DIRECTORY)
            install(FILES "${_f}" DESTINATION "petakit5d_fftw/licenses/${_sub}")
        endforeach()
    endif()
endif()
