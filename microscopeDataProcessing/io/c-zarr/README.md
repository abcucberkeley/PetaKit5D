# c-zarr
An efficient parallel Zarr reader/writer that utilizes c-blosc/c-blosc2 and OpenMP.

## Quick Start Guide (MATLAB)

### Prerequisites
1. (Windows Only) All neccessary libraries are included for the windows version.
2. (Linux Only) You will need c-blosc, c-blosc2, and cJSON which can all be obtained and built with the included buildc-zarr.sh script.


### Download and Install
1. Download the latest release for your OS from here (windows/linux.zip): https://github.com/abcucberkeley/c-zarr/releases
2. Unzip the folder
3. (Windows Only) You can now put the folders wherever you'd like and add them to your path if needed. Just keep the mex files with their associated dll files so the function can always run.
4. (Linux Only) chmod +x buildc-zarr.sh to make it executable. You can now run the buildc-zarr.sh script and follow the instructions. Now when you open matlab, just add the mex files to your path when needed

### Usage
Note: For filepath separators, on Mac/Linux you can use / and on Windows you can use \

#### createZarrFile - Create a custom .zarray metadata file
% Note the created file is probably hidden by default on your system

createZarrFile('path/to/file.zarr');

#### parallelReadZarr - Read a Zarr image into an array
im = parallelReadZarr('path/to/file.zarr');

#### parallelWriteTiff - Write an array out as a Zarr image
im = rand(100,100,100);

% The third input can always be 1 to use a uuid for the written blocks

% The fourth input is the size of the blocks

parallelWriteZarr('path/to/file.zarr',im, 1, [256,256,256]);
