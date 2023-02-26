# c-zarr
An efficient parallel Zarr reader/writer that utilizes c-blosc/c-blosc2 and OpenMP.

## Quick Start Guide (MATLAB)

### Prerequisites
1. All neccessary libraries are included for the Linux and Windows version.

### Download and Install
1. Download the latest release for your OS from here (windows.zip/linux.tar.gz): https://github.com/abcucberkeley/c-zarr/releases
2. Unzip the folder
3. You can now put the folders wherever you'd like and add them to your path if needed. Just keep the mex files with their associated dll files so the mex function can always run.

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
