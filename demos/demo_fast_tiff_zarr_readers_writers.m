%% demo to use our fast tiff/zarr readers and writers
% also include comparison with conventional ones
% 
% Note: the acceleration over convetional methods depends on the number of
% CPU cores in your system. In the paper, we benchmarked in a 24-core
% system. If your system has more cores, you may see even more accelerations. 
% However, if your system has fewer cores, the acceleration may be not as much 
% as we shown in the paper. 

clear, clc;
setup();

%% Step 1: get our demo data from zenodo/Dropbox (skip this step if the data is already downloaded)
% download the example dataset from zenodo (https://doi.org/10.5281/zenodo.10471978) manually, 
% or use the code below to download the data from Dropbox
if ispc
    destPath = fullfile(getenv('USERPROFILE'), 'Downloads');
    destPath = strrep(destPath, '\', '/');
else
    destPath = '~/Downloads/';
end
demo_data_downloader(destPath);

dataPath = [destPath, '/LLSM5DTools_demo_cell_image_dataset/'];


%% read a tiff file

fsn = 'Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_003y_000z_0000t';
fn = [dataPath, fsn, '.tif'];

% conventional tiff reader
fprintf('Conventional Tiff reader (libtiff): ');
tic
im_1 = readtiff_matlab(fn);
toc
clear im_1;

% fast Tiff reader
fprintf('Cpp-Tiff reader: ');
tic
im = readtiff(fn);
toc
fprintf('\n');

%% replicate the data to 10000 frames
% 
nframe = 5000;

im_rep = repmat(im, 1, 1, round(nframe / 500));
disp(size(im_rep))

outPath = [dataPath, 'replicated/'];
mkdir(outPath);


%% write the replicate data as Tiff format to disk

% write to disk as a Tiff file

fnout = sprintf('%s%s_frame_number_%d_conventional.tif', outPath, fsn, nframe);
% conventional Tiff writer
fprintf('Conventional Tiff writer (libtiff): ');
tic
writetiff(im_rep, fnout, mode='libtiff');
toc

fnout = sprintf('%s%s_frame_number_%d.tif', outPath, fsn, nframe);
% fast Tiff writer
fprintf('Cpp-Tiff writer: ');
tic
writetiff(im_rep, fnout);
toc
fprintf('\n');


%% read the replicated data

clear im_rep;

% conventional Tiff reader
fprintf('Conventional Tiff reader (libtiff) for %d frames: ', nframe);
tic
im_rep = readtiff_matlab(fnout);
toc
clear im_rep;

% fast Tiff reader
fprintf('Cpp-Tiff reader for %d frames: ', nframe);
tic
im_rep = readtiff(fnout);
toc
fprintf('\n');


%% write the replicate data as Zarr format to disk

zarrFnout = sprintf('%s%s_frame_number_%d.zarr', outPath, fsn, nframe);

% create Zarr file
dataSize = size(im_rep);
% chunk size for Zarr
blockSize = [256, 256, 256];
createzarr(zarrFnout, dataSize=dataSize, blockSize=blockSize);

% write Zarr with our faster writer
fprintf('Cpp-zarr writer for %d frames: ', nframe);
tic
writezarr(im_rep, zarrFnout);
toc
clear im_rep;

%% read Zarr with our faster writer

fprintf('Cpp-zarr reader for %d frames: ', nframe);
tic
im_rep = readzarr(zarrFnout);
toc


%% read/write Zarr with conventional zarr library (MATLAB interface of Zarr)
% note: to make it work, conda in python need to be installed, and Zarr needs 
% to be installed in an enviroment. Then, the environment need to be
% activated and Matlab is launched in that enviroment

% setup environment to load python environment for matlab
setup([], true);

zarrFnout_1 = sprintf('%s%s_frame_number_%d_conventional.zarr', outPath, fsn, nframe);

% create zarr file
dtype = class(im_rep);
dataSize = size(im_rep);
blockSize = [256, 256, 256];

init_val = zeros(1, dtype);
bim = blockedImage(zarrFnout_1, dataSize, blockSize, init_val, "Adapter", ZarrAdapter, 'Mode', 'w');             
bim.Adapter.close();

% write zarr with conventional Zarr writer
fprintf('MATLAB interface of Zarr writer for %d frames: ', nframe);
tic
bim = blockedImage(zarrFnout_1, "Adapter", ZarrAdapter);
bim.Adapter.setRegion(ones(1, numel(bim.Size)), bim.Size, im_rep)
toc


%% read zarr with conventional method

clear im_rep;
fprintf('MATLAB interface of Zarr reader for %d frames: ', nframe);
tic
bim = blockedImage(zarrFnout_1, "Adapter", ZarrAdapter);
im_rep = bim.Adapter.getIORegion(ones(1, numel(bim.Size)), bim.Size);
toc




















