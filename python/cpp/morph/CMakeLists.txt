# petakit5d._ext.morph: binary morphology as MATLAB's morphop applies it, the exact Euclidean
# distance transform (bwdist), connected components in MATLAB's orders, convex hull images
# (bwconvhull) and median filters (medfilt2, medfilt3)
pk5d_add_extension(morph module.cpp binary.cpp edt.cpp hull.cpp label.cpp median.cpp)
