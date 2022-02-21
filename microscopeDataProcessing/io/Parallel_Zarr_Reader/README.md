# Parallel_Zarr_Reader
Note: cJSON and c-blosc2 are required
# MATLAB USAGE

The mex files can be called directly from MATLAB once added to the MATLAB path.

Windows mex extentsion - .mexw64

Linux mex extentsion - .mexa64

Mac mex extentsion - .mexmaci64

EXAMPLE USE (Returns a 3D array in y,x,z format):

myZarr = parallelReadZarr('C:\Users\Example\Desktop\test.zarr');

USER DEFINED RANGE EXAMPLE USE (Second argument is a 3D bounding box [startY, startX, startZ, endY, endX, endZ]):

myZarr = parallelReadZarr('C:\Users\Example\Desktop\test.zarr',[1, 1, 1, 1000, 1000, 1000]);
