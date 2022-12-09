# c-tiff
An efficient parallel Tiff reader/writer that utilizes LibTIFF and OpenMP.

## Quick Start Guide (MATLAB)

### Prerequisites
1. None! The parallel reader and writer mex files will work with most recent version of Matlab.

### Download and Install
1. Download the latest release for your OS from here (windows/linux/mac.zip): https://github.com/abcucberkeley/c-tiff/releases
2. Unzip the folder
3. You can now put the mex files wherever you'd like and add them to your path if needed

### Usage
Note: For filepath separators, on Mac/Linux you can use / and on Windows you can use \

#### parallelReadTiff - Read a Tiff image into an array
im = parallelReadTiff('path/to/file.tif');

#### getImageSize_mex - Get the dimensions of a Tiff image
size = getImageSize_mex('path/to/file.tif');

#### parallelWriteTiff - Write an array out as a Tiff image
im = rand(100,100,100);

parallelWriteTiff('path/to/file.tif',im);
