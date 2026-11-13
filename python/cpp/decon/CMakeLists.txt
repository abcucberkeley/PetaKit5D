# petakit5d._ext.decon: the deconvolution kernels moved from the mex files of
# microscopeDataProcessing/decon/mex: the half-spectrum Richardson-Lucy steps of decon_lucy_rfft_mex (in
# MATLAB's layout on an FFT backend's table, and on the port's layout with the native half spectra), the
# acceleration step and fixed-order sums of decon_lucy_accel_mex, the mask erosion of
# decon_mask_edge_erosion_mex and the window bounds, components and sort of omw_backprojector_mex. The FFTs
# come from cpp/fft's driver library (no FFT library header is included here).
pk5d_add_extension(decon
  module.cpp
  transpose.cpp
  rfft_matlab.cpp
  rfft_native.cpp
  accel.cpp
  erosion.cpp
  omw.cpp
  LIBRARIES pk5d_fft_driver)
