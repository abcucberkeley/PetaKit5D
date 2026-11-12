# FFTW 3.3.x (GPL-2.0-or-later), single and double precision, as two shared libraries with FFTW's API
# exported, built with the codelet sets of MATLAB R2024b's libmwfftw3f / libmwfftw3 ("fftw-3.3.8-sse2-avx":
# SSE2 and AVX, no FMA): the port's optional exact FFT backend, loaded at run time through the backend table
# (cpp/fft/backend.cpp). The release tarball is fetched by URL and SHA-256 (FETCHCONTENT_SOURCE_DIR_PK5D_FFTW
# points to an unpacked copy instead). The sources are compiled from FFTW's own file lists; as FFTW's configure
# does, the SIMD flags go only on the SIMD codelets, and FFTW picks them at run time from the CPU. FFTW's
# threads and Fortran wrappers are left out.
#
# Defines the targets pk5d_fftw3f and pk5d_fftw3 (output in PK5D_FFTW_OUTPUT_DIR).

set(PK5D_FFTW_VERSION "3.3.8" CACHE STRING "FFTW release: 3.3.8 (MATLAB R2024b's, default) or 3.3.10")
set_property(CACHE PK5D_FFTW_VERSION PROPERTY STRINGS 3.3.8 3.3.10)
set(_pk5d_fftw_sha256_3.3.8 6113262f6e92c5bd474f2875fa1b01054c4ad5040f6b0da7c03c98821d9ae303)
set(_pk5d_fftw_sha256_3.3.10 56c932549852cddcfafdab3820b0200c7742675be92179e59e6215b340e26467)
if(NOT DEFINED _pk5d_fftw_sha256_${PK5D_FFTW_VERSION})
    message(FATAL_ERROR "PK5D_FFTW_VERSION must be 3.3.8 or 3.3.10")
endif()
set(PK5D_FFTW_URL "https://www.fftw.org/fftw-${PK5D_FFTW_VERSION}.tar.gz" CACHE STRING "FFTW tarball (any mirror)")

include(FetchContent)
FetchContent_Declare(pk5d_fftw
    URL ${PK5D_FFTW_URL}
    URL_HASH SHA256=${_pk5d_fftw_sha256_${PK5D_FFTW_VERSION}}
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    # sources only: FFTW's own CMakeLists.txt is not used
    SOURCE_SUBDIR pk5d-no-cmake)
FetchContent_MakeAvailable(pk5d_fftw)
set(_fs ${pk5d_fftw_SOURCE_DIR})

# --- configuration, with FFTW's own checks ---
include(CheckIncludeFile)
include(CheckSymbolExists)
include(CheckTypeSize)
include(CheckPrototypeDefinition)
include(CheckCCompilerFlag)
foreach(_h alloca.h dlfcn.h inttypes.h limits.h malloc.h memory.h stddef.h stdint.h stdlib.h string.h strings.h
        sys/types.h sys/time.h sys/stat.h time.h unistd.h)
    string(TOUPPER "HAVE_${_h}" _v)
    string(REGEX REPLACE "[./]" "_" _v ${_v})
    check_include_file(${_h} ${_v})
endforeach()
if(HAVE_TIME_H AND HAVE_SYS_TIME_H)
    set(TIME_WITH_SYS_TIME TRUE)
endif()
check_prototype_definition(drand48 "double drand48 (void)" "0" stdlib.h HAVE_DECL_DRAND48)
check_prototype_definition(srand48 "void srand48(long int seedval)" "0" stdlib.h HAVE_DECL_SRAND48)
check_prototype_definition(cosl "long double cosl( long double arg )" "0" math.h HAVE_DECL_COSL)
check_prototype_definition(sinl "long double sinl( long double arg )" "0" math.h HAVE_DECL_SINL)
check_prototype_definition(memalign "void *memalign(size_t alignment, size_t size)" "0" malloc.h HAVE_DECL_MEMALIGN)
check_prototype_definition(posix_memalign "int posix_memalign(void **memptr, size_t alignment, size_t size)" "0" stdlib.h
    HAVE_DECL_POSIX_MEMALIGN)
check_symbol_exists(clock_gettime time.h HAVE_CLOCK_GETTIME)
check_symbol_exists(gettimeofday sys/time.h HAVE_GETTIMEOFDAY)
check_symbol_exists(getpagesize unistd.h HAVE_GETPAGESIZE)
check_symbol_exists(drand48 stdlib.h HAVE_DRAND48)
check_symbol_exists(srand48 stdlib.h HAVE_SRAND48)
check_symbol_exists(memalign malloc.h HAVE_MEMALIGN)
check_symbol_exists(posix_memalign stdlib.h HAVE_POSIX_MEMALIGN)
check_symbol_exists(mach_absolute_time mach/mach_time.h HAVE_MACH_ABSOLUTE_TIME)
check_symbol_exists(alloca alloca.h HAVE_ALLOCA)
if(NOT HAVE_ALLOCA)
    unset(HAVE_ALLOCA CACHE)
    check_symbol_exists(alloca malloc.h HAVE_ALLOCA)
endif()
check_symbol_exists(isnan math.h HAVE_ISNAN)
check_symbol_exists(snprintf stdio.h HAVE_SNPRINTF)
check_symbol_exists(sysctl unistd.h HAVE_SYSCTL)
if(UNIX)
    set(CMAKE_REQUIRED_LIBRARIES m)
endif()
check_symbol_exists(cosl math.h HAVE_COSL)
check_symbol_exists(sinl math.h HAVE_SINL)
unset(CMAKE_REQUIRED_LIBRARIES)
foreach(_t float double int long "long long" "unsigned int" "unsigned long" "unsigned long long" size_t ptrdiff_t)
    string(TOUPPER "SIZEOF_${_t}" _v)
    string(REPLACE " " "_" _v ${_v})
    check_type_size("${_t}" ${_v})
endforeach()
if(UNIX)
    set(HAVE_LIBM TRUE)
endif()

# SSE2 and AVX on x86-64 only (the codelet sets of MATLAB's build); other systems get the scalar codelets
set(_pk5d_x86 FALSE)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64|x64|i[3-6]86)$")
    set(_pk5d_x86 TRUE)
endif()
set(_pk5d_sse2_flag "")
set(_pk5d_avx_flag "")
if(_pk5d_x86)
    if(MSVC)
        check_c_compiler_flag("/arch:AVX" PK5D_FFTW_HAVE_ARCH_AVX)
        set(_pk5d_avx_flag "/arch:AVX")
    else()
        check_c_compiler_flag("-msse2" PK5D_FFTW_HAVE_MSSE2)
        check_c_compiler_flag("-mavx" PK5D_FFTW_HAVE_MAVX)
        set(_pk5d_sse2_flag "-msse2")
        set(_pk5d_avx_flag "-mavx")
    endif()
    set(HAVE_SSE2 TRUE)
    set(HAVE_AVX TRUE)
endif()
set(FFTW_VERSION ${PK5D_FFTW_VERSION})
set(DISABLE_FORTRAN TRUE)

# FFTW's config.h template, with a compiler name instead of the compiler's path in fftw_cc
file(READ ${_fs}/cmake.config.h.in _pk5d_cfg)
string(REPLACE "@CMAKE_C_COMPILER@" "${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}" _pk5d_cfg "${_pk5d_cfg}")
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/fftw_config.h.in "${_pk5d_cfg}")
set(FFTW_SINGLE TRUE)
configure_file(${CMAKE_CURRENT_BINARY_DIR}/fftw_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/fftw-single/config.h @ONLY)
unset(FFTW_SINGLE)
configure_file(${CMAKE_CURRENT_BINARY_DIR}/fftw_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/fftw-double/config.h @ONLY)

# --- sources: FFTW's lists without threads, Fortran wrappers and benchmarks ---
file(GLOB _pk5d_fftw_common
    ${_fs}/api/*.c ${_fs}/dft/*.c ${_fs}/dft/scalar/*.c ${_fs}/dft/scalar/codelets/*.c ${_fs}/dft/simd/*.c
    ${_fs}/kernel/*.c ${_fs}/rdft/*.c ${_fs}/rdft/scalar/*.c ${_fs}/rdft/scalar/r2cb/*.c ${_fs}/rdft/scalar/r2cf/*.c
    ${_fs}/rdft/scalar/r2r/*.c ${_fs}/rdft/simd/*.c ${_fs}/reodft/*.c ${_fs}/simd-support/*.c)
list(REMOVE_ITEM _pk5d_fftw_common ${_fs}/api/f77api.c)
set(_pk5d_fftw_sse2 "")
set(_pk5d_fftw_avx "")
if(_pk5d_x86)
    file(GLOB _pk5d_fftw_sse2 ${_fs}/dft/simd/sse2/*.c ${_fs}/rdft/simd/sse2/*.c)
    file(GLOB _pk5d_fftw_avx ${_fs}/dft/simd/avx/*.c ${_fs}/rdft/simd/avx/*.c)
    if(_pk5d_sse2_flag)
        set_source_files_properties(${_pk5d_fftw_sse2} PROPERTIES COMPILE_OPTIONS "${_pk5d_sse2_flag}")
    endif()
    set_source_files_properties(${_pk5d_fftw_avx} PROPERTIES COMPILE_OPTIONS "${_pk5d_avx_flag}")
endif()

foreach(_prec single double)
    if(_prec STREQUAL "single")
        set(_t pk5d_fftw3f)
    else()
        set(_t pk5d_fftw3)
    endif()
    add_library(${_t} SHARED ${_pk5d_fftw_common} ${_pk5d_fftw_sse2} ${_pk5d_fftw_avx})
    target_include_directories(${_t} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/fftw-${_prec} ${_fs})
    set_target_properties(${_t} PROPERTIES POSITION_INDEPENDENT_CODE ON C_VISIBILITY_PRESET default
        LIBRARY_OUTPUT_DIRECTORY ${PK5D_FFTW_OUTPUT_DIR} RUNTIME_OUTPUT_DIRECTORY ${PK5D_FFTW_OUTPUT_DIR})
    if(MSVC)
        # no fused multiply-adds, IEEE semantics; FFTW's API exported from the DLL
        target_compile_options(${_t} PRIVATE /fp:precise /bigobj)
        target_compile_definitions(${_t} PRIVATE _CRT_SECURE_NO_WARNINGS FFTW_DLL)
    else()
        target_compile_options(${_t} PRIVATE -ffp-contract=off -fno-fast-math -w)
    endif()
    if(UNIX)
        target_link_libraries(${_t} PRIVATE m)
    endif()
endforeach()
