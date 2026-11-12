# petakit5d._ext.fft: MATLAB R2024b's fft, ifft, fftn and ifftn, bit for bit, on our own FFTW build.
# Library targets for other groups (LIBRARIES of pk5d_add_extension):
#   pk5d_fftw        FFTW 3.3.x single and double precision, static, symbols hidden (fftw.cmake)
#   pk5d_fft_driver  MATLAB's FFT driver on it (driver.h), with the test-only switch to another FFTW
include(${CMAKE_CURRENT_LIST_DIR}/fftw.cmake)

add_library(pk5d_fftw INTERFACE)
target_link_libraries(pk5d_fftw INTERFACE pk5d_fftw3f pk5d_fftw3)

add_library(pk5d_fft_driver STATIC backend.cpp driver.cpp)
target_include_directories(pk5d_fft_driver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(pk5d_fft_driver PUBLIC pk5d_fftw PRIVATE pk5d_flags pk5d_openmp)
set_target_properties(pk5d_fft_driver PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_VISIBILITY_PRESET hidden)
if(NOT WIN32)
    target_link_libraries(pk5d_fft_driver PUBLIC ${CMAKE_DL_LIBS})
endif()

pk5d_add_extension(fft module.cpp LIBRARIES pk5d_fft_driver)
