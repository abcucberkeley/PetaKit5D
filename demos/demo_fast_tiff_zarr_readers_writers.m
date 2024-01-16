%% demo to use our fast tiff/zarr readers and writers
% also include comparison with conventional ones
% 
% Note: the acceleration over conventional methods depends on the number of
% CPU cores in your system. In the paper, we benchmarked in a 24-core
% system. If your system has more cores, you may see even more accelerations. 
% However, if your system has fewer cores, the acceleration may be not as much 
% as we showed in the paper. 

clear, clc;

fprintf('Fast Tiff/Zarr readers and writers demo...\n\n');

% move to the LLSM5DTools root directory
curPath = pwd;
if ~endsWith(curPath, 'LLSM5DTools')
    mfilePath = mfilename('fullpath');
    if contains(mfilePath,'LiveEditorEvaluationHelper')
        mfilePath = matlab.desktop.editor.getActiveFilename;
    end
    
    mPath = fileparts(mfilePath);
    if endsWith(mPath, 'demos')
        cd(mPath);
        cd('..')
    end
end

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
fprintf('Conventional Tiff reader (libtiff) for 500 frames: ');
tic
im_1 = readtiff_matlab(fn);
toc
clear im_1;

% fast Tiff reader
% if im already exist, clear it to avoid the overhead to overwrite the variable
if exist('im', 'var')
    clear im;
end
fprintf('Cpp-Tiff reader for 500 frames: ');
tic
im = readtiff(fn);
toc
fprintf('\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% benchmark readers and writers on larger data with more frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% replicate the data to 10000 frames 

% Note: please make sure your system has at least 64 GB available RAM for
% 10k frames; otherwise, you may need to reduce the number of frames. 

nframe = 10000;

im_rep = repmat(im, 1, 1, ceil(nframe / 500));
if size(im_rep, 3) > nframe
    im_rep = crop3d_mex(im_rep, [1, 1, 1, size(im_rep, 1), size(im_rep, 2), nframe]);
end
disp(size(im_rep));

outPath = [dataPath, 'replicated/'];
mkdir(outPath);


%% write the replicated data as Tiff format to disk

% write to disk as a Tiff file

fnout = sprintf('%s%s_frame_number_%d_conventional.tif', outPath, fsn, nframe);
% conventional Tiff writer
fprintf('Conventional Tiff writer (libtiff) for %d frames (takes ~3 min for 16-core CPUs): ', nframe);
tic
writetiff(im_rep, fnout, mode='libtiff');
toc

fnout = sprintf('%s%s_frame_number_%d.tif', outPath, fsn, nframe);
% fast Tiff writer
fprintf('Cpp-Tiff writer for %d frames: ', nframe);
tic
writetiff(im_rep, fnout);
toc
fprintf('\n');


%% read the replicated data

clear im_rep;

% conventional Tiff reader
fprintf('Conventional Tiff reader (libtiff) for %d frames (takes ~1.5 min for 16-core CPUs): ', nframe);
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


%% write the replicated data as Zarr format to disk

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

%% read Zarr with our fast reader

fprintf('Cpp-zarr reader for %d frames: ', nframe);
tic
im_rep = readzarr(zarrFnout);
toc


%% read/write Zarr with conventional zarr library (MATLAB interface of Zarr)
% Note: to make it work, anaconda in python need to be installed, and
% zarr-python and dask need to be installed in an enviroment. Then, the environment 
% needs to be activated and Matlab is launched from that enviroment.

% this section and the following one is disabled by default. If you
% installed Python enviroment as described above, please comment out the
% code belew. 
fprintf(['\nThe conventional Zarr reader and writer requies a Python enviroment as desribed in the comment above.\n', ...
         'If you have installed the Python enviroment and launched Matlab from the environment,please comment \n', ...
         'out the line below to allow the running of the convetional Zarr reader and writer.\n']);
return;

% setup environment to load the python environment for matlab
% please provide your python executable path, i.e.,
%   ~/anaconda3/envs/your_env/bin/python for linux and MacOS
%   C:\Users\UserName\anaconda3\envs\your_env\python.EXE for Windows
pythonPath = '/path/to/your/python';
setup([], true, pythonPath);

zarrFnout_1 = sprintf('%s%s_frame_number_%d_conventional.zarr', outPath, fsn, nframe);

% create zarr file
dtype = class(im_rep);
dataSize = size(im_rep);
blockSize = [256, 256, 256];

init_val = zeros(1, dtype);
if exist(zarrFnout_1, 'dir')
    rmdir(zarrFnout_1, "s");
end
bim = blockedImage(zarrFnout_1, dataSize, blockSize, init_val, "Adapter", ZarrAdapter, 'Mode', 'w');             
bim.Adapter.close();

% write zarr with conventional Zarr writer
fprintf('MATLAB interface of Zarr writer for %d frames: ', nframe);
tic
bim = blockedImage(zarrFnout_1, "Adapter", ZarrAdapter);
bim.Adapter.setRegion(ones(1, numel(bim.Size)), bim.Size, im_rep)
toc


%% read zarr with the conventional method
% also requires the Python environment as metioned above

clear im_rep;
fprintf('MATLAB interface of Zarr reader for %d frames: ', nframe);
tic
bim = blockedImage(zarrFnout_1, "Adapter", ZarrAdapter);
im_rep = bim.Adapter.getIORegion(ones(1, numel(bim.Size)), bim.Size);
toc




