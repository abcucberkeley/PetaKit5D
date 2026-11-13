# petakit5d._ext.crop: the box copies of crop3d_mex, crop4d_mex, indexing3d_mex, indexing4d_mex
# and indexing4d_crop_mex (microscopeDataProcessing/crop/mex); the per-region thread binding
# comes from the element-wise primitives (core/ops/ops_bind.hpp, header-only)
pk5d_add_extension(crop module.cpp INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../core/ops")
