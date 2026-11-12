# petakit5d._ext.deskew: the deskew/rotate warp and the skewed space interpolation kernels
# (the three mex files of microscopeDataProcessing/deskew_rotate/mex). -mavx2 on x86-64, as the
# mex are built: the line sums use the mex's AVX2 code there and the same order elsewhere.
pk5d_add_extension(deskew module.cpp AVX2)

# petakit5d._ext.deskew_scalar: the same sources without AVX2 (the code other processors
# run), so that x86 builds and CI test it against the AVX2 module bit for bit
if(PK5D_X86_64)
  option(PK5D_DESKEW_SCALAR "Build petakit5d._ext.deskew_scalar (the deskew kernels without AVX2)" ON)
  if(PK5D_DESKEW_SCALAR)
    pk5d_add_extension(deskew_scalar module.cpp)
  endif()
endif()
