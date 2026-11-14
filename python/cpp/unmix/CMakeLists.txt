# petakit5d._ext.unmix: the weighted channel sum of unmix_channels_sum_mex
# (microscopeDataProcessing/unmix/mex); the per-region thread binding comes from the
# element-wise primitives (core/ops/ops_bind.hpp, header-only)
pk5d_add_extension(unmix module.cpp INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../core/ops")
