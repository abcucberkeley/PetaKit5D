# petakit5d._ext.flatfield: the flat-field correction of flat_field_correction_3d_mex
# (microscopeDataProcessing/flatfield_correction/mex); the per-region thread binding comes
# from the element-wise primitives (core/ops/ops_bind.hpp, header-only)
pk5d_add_extension(flatfield module.cpp INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../core/ops")
