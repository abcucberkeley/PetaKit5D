# petakit5d._ext.filters: the arithmetic of MATLAB's linear filters (N04): tap programs along one
# axis with MATLAB's padding (the passes of imfilter and imgaussfilt) and N-D correlations in a
# given tap order (convn, conv2), with fused multiply-adds where MATLAB's code has them
pk5d_add_extension(filters module.cpp)
