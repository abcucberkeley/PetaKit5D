# petakit5d._ext.xcorr: the registration's normalized cross-correlation in a window of its full output
# (normxcorr3_fast and normxcorr3_window_mex), on the FFT driver and backends of the fft group
pk5d_add_extension(xcorr module.cpp window.cpp ncc.cpp LIBRARIES pk5d_fft_driver)
