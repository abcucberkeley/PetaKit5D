# The CPU FFT: MATLAB's fft driver and the half-spectrum transforms (BSD-3), with pocketfft (BSD-3,
# third_party/pocketfft) built in as the default ("permissive") backend. FFTW 3.3.8 with MATLAB's codelets
# (GPL), the "exact" backend, is a separate add-on of two shared libraries loaded at run time (exact/); the
# default build neither builds nor links it, and no file of it includes fftw3.h.
add_library(pk5d_fft_driver STATIC backend.cpp driver.cpp pocketfft_api.cpp halfspec.cpp)
target_include_directories(pk5d_fft_driver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
                           PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/third_party/pocketfft)
target_link_libraries(pk5d_fft_driver PRIVATE pk5d_flags pk5d_openmp)
set_target_properties(pk5d_fft_driver PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_VISIBILITY_PRESET hidden)
find_package(Threads REQUIRED)
target_link_libraries(pk5d_fft_driver PUBLIC Threads::Threads)
if(NOT WIN32)
    target_link_libraries(pk5d_fft_driver PUBLIC ${CMAKE_DL_LIBS})
endif()
if(MSVC AND NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # pocketfft's vector types need GCC's or Clang's vector extensions: with MSVC every line runs alone (the
    # same bits, slower); build with clang-cl for the vectors (docs/decisions/0006-fft-backends.md)
    message(STATUS "petakit5d fft: MSVC builds pocketfft without vectors; clang-cl builds it with them")
endif()

pk5d_add_extension(fft module.cpp LIBRARIES pk5d_fft_driver)

# The exact add-on, for development and the oracle tests: FFTW 3.3.8 (fetched by URL and SHA-256) as two
# shared libraries in <build>/fftw_exact, never installed with the package. Its own package: exact/pyproject.toml.
option(PK5D_FFT_EXACT "Also build the exact FFT add-on (FFTW 3.3.8, GPL) into <build>/fftw_exact" OFF)
if(PK5D_FFT_EXACT)
    if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/exact/CMakeLists.txt)
        # the package's source distribution leaves the add-on out
        message(FATAL_ERROR "PK5D_FFT_EXACT: the add-on's recipe (cpp/fft/exact) is not in this source tree; "
                            "build petakit5d-fftw from the repository")
    endif()
    add_subdirectory(exact)
endif()
