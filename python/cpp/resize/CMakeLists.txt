# petakit5d._ext.resize: the passes of MATLAB's imresize and imresize3 (N01): weighted sums along
# one axis in double with the rounding of MATLAB's resize, index copies for one-hot weights, and
# the cubic kernel (fdlibm's pow); baseline and AVX2 code paths chosen at run time
pk5d_add_extension(resize module.cpp)
