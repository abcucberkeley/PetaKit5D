# petakit5d._ext.stitch: the stitching kernels (feather blending, distance map resize, tile
# blending, integral image, the major-tile check, NaN replacement), moved from the mex files of
# microscopeDataProcessing/stitch/mex and pointDetection/mex_functions/replace_nan_with_value.c
pk5d_add_extension(stitch
  module.cpp
  feather_blending.cpp
  distance_map_resize.cpp
  tile_blending.cpp
  integral_image.cpp
  any_4th_dim.cpp
  replace_nan.cpp)
