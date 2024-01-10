% demo to run RL deconvolution with conventional and OMW methods


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


%% Step 2: test the parameter for OMW backward projector

% Cam A
psfFn = [dataPath, 'PSF/560_c.tif'];
% OTF thresholding parameter
OTFCumThresh = 0.9;
% true if the PSF is in skew space
skewed = true;
XR_visualize_OTF_mask_segmentation(psfFn, OTFCumThresh, skewed);


%% 

% Cam B
psfFn = [dataPath, 'PSF/488_2_c.tif'];
% OTF thresholding parameter
OTFCumThresh = 0.9;
% true if the PSF is in skew space
skewed = true;
XR_visualize_OTF_mask_segmentation(psfFn, OTFCumThresh, skewed);


%% Step 3: OMW deconvolution 
%% Step 3.1: set parameters 
% add the software to the path
setup([]);

% root path
rt = dataPath;
% data path for data to be deconvolved, also support for multiple data folders
dataPaths = {[rt, '/']};

% xy pixel size in um
xyPixelSize = 0.108;
% z step size
dz = 0.3;
% scan direction
Reverse = true;
% psf z step size (we assume xyPixelSize also apply to psf)
dzPSF = 0.3;

% if true, check whether image is flipped in z using the setting files
parseSettingFile = false;

% channel patterns for the channels, the channel patterns should map the
% order of PSF filenames.
ChannelPatterns = {'CamB_ch0', ...
                   };  

% psf path
psf_rt = rt;            
PSFFullpaths = {
                [psf_rt, 'PSF/488_2_c.tif'], ...
                };            

% RL method
RLmethod = 'omw';
% wiener filter parameter
wienerAlpha = 0.005;
% OTF thresholding parameter
OTFCumThresh = 0.9;
% true if the PSF is in skew space
skewed = true;
% deconvolution result path string (within dataPath)
deconPathstr = 'matlab_decon_omw';

% background to subtract
Background = 100;
% number of iterations
DeconIter = 2;
% decon to 80 iterations (not use the criteria for early stop)
fixIter = true;
% erode the edge after decon for number of pixels.
EdgeErosion = 8;
% save as 16bit, if false, save to single
Save16bit = true;
% use zarr file as input, if false, use tiff as input
zarrFile = false;
% number of cpu cores
cpusPerTask = 24;
% use cluster computing for different images
parseCluster = false;
% set it to true for large files that cannot be fitted to RAM/GPU, it will
% split the data to chunks for deconvolution
largeFile = false;
% use GPU for deconvolution
GPUJob = true;
% if true, save intermediate results every 5 iterations.
debug = false;

%% Step 3.2: run the deconvolution with given parameters. 
% the results will be saved in matlab_decon under the dataPaths. 
% the next step is deskew/rotate (if in skewed space for x-stage scan) or 
% rotate (if objective scan) or other processings. 

% result folder:
% {destPath}/LLSM5DTools_demo_cell_image_dataset/matlab_decon_omw/

XR_decon_data_wrapper(dataPaths, 'deconPathstr', deconPathstr, 'xyPixelSize', xyPixelSize, ...
    'dz', dz, 'Reverse', Reverse, 'ChannelPatterns', ChannelPatterns, 'PSFFullpaths', PSFFullpaths, ...
    'dzPSF', dzPSF, 'parseSettingFile', parseSettingFile, 'RLmethod', RLmethod, ...
    'wienerAlpha', wienerAlpha, 'OTFCumThresh', OTFCumThresh, 'skewed', skewed, ...
    'Background', Background, 'CPPdecon', false, 'CudaDecon', false, 'DeconIter', DeconIter, ...
    'fixIter', fixIter, 'EdgeErosion', EdgeErosion, 'Save16bit', Save16bit, ...
    'zarrFile', zarrFile, 'parseCluster', parseCluster, 'largeFile', largeFile, ...
    'GPUJob', GPUJob, 'debug', debug, 'cpusPerTask', cpusPerTask);


%% Step 4: Conventional RL deconvolution 
%% Step 4.1: set parameters 
% add the software to the path
setup([]);

% root path
rt = dataPath;
% data path for data to be deconvolved, also support for multiple data folders
dataPaths = {[rt, '/']};

% xy pixel size in um
xyPixelSize = 0.108;
% z step size
dz = 0.3;
% scan direction
Reverse = true;
% psf z step size (we assume xyPixelSize also apply to psf)
dzPSF = 0.3;

% if true, check whether image is flipped in z using the setting files
parseSettingFile = false;

% channel patterns for the channels, the channel patterns should map the
% order of PSF filenames.
ChannelPatterns = {'CamB_ch0', ...
                   };  

% psf path
psf_rt = rt;            
PSFFullpaths = {
                [psf_rt, 'PSF/488_2_c.tif'], ...
                };            

% RL method
RLmethod = 'simplified';
% deconvolution result path string (within dataPath)
deconPathstr = 'matlab_decon_conventional';

% background to subtract
Background = 100;
% number of iterations
DeconIter = 30;
% decon to 80 iterations (not use the criteria for early stop)
fixIter = true;
% erode the edge after decon for number of pixels.
EdgeErosion = 8;
% save as 16bit, if false, save to single
Save16bit = true;
% use zarr file as input, if false, use tiff as input
zarrFile = false;
% number of cpu cores
cpusPerTask = 24;
% use cluster computing for different images
parseCluster = false;
% set it to true for large files that cannot be fitted to RAM/GPU, it will
% split the data to chunks for deconvolution
largeFile = false;
% use GPU for deconvolution
GPUJob = true;
% if true, save intermediate results every 5 iterations.
debug = false;

%% Step 4.2: run the deconvolution with given parameters. 
% the results will be saved in matlab_decon under the dataPaths. 
% the next step is deskew/rotate (if in skewed space for x-stage scan) or 
% rotate (if objective scan) or other processings. 

% result folder:
% {destPath}/LLSM5DTools_demo_cell_image_dataset/matlab_decon_conventional/

XR_decon_data_wrapper(dataPaths, 'deconPathstr', deconPathstr, 'xyPixelSize', xyPixelSize, ...
    'dz', dz, 'Reverse', Reverse, 'ChannelPatterns', ChannelPatterns, 'PSFFullpaths', PSFFullpaths, ...
    'dzPSF', dzPSF, 'parseSettingFile', parseSettingFile, 'RLmethod', RLmethod, ...
    'wienerAlpha', wienerAlpha, 'OTFCumThresh', OTFCumThresh, 'skewed', skewed, ...
    'Background', Background, 'CPPdecon', false, 'CudaDecon', false, 'DeconIter', DeconIter, ...
    'fixIter', fixIter, 'EdgeErosion', EdgeErosion, 'Save16bit', Save16bit, ...
    'zarrFile', zarrFile, 'parseCluster', parseCluster, 'largeFile', largeFile, ...
    'GPUJob', GPUJob, 'debug', debug, 'cpusPerTask', cpusPerTask);


%% Step 5: deskew/rotate the deconvoluved results

% result folders:
% {destPath}/LLSM5DTools_demo_cell_image_dataset/DSR/
% {destPath}/LLSM5DTools_demo_cell_image_dataset/matlab_decon_omw/DSR/
% {destPath}/LLSM5DTools_demo_cell_image_dataset/matlab_decon_conventional/DSR/

dataPath_exps = {
                 [dataPath], ...
                 [dataPath,'matlab_decon_omw/'], ...
                 [dataPath,'matlab_decon_conventional/'], ...
                };

% xy pixel size in um
xyPixelSize = 0.108;
% z step size
dz = 0.3;
% scan direction
Reverse = true;
% channel patterns for the channels to process
ChannelPatterns = {'CamB_ch0', ...
                   };  

% true if using large-scale process method
largeFile = false;
% true if input is in Zarr format
zarrFile = false;
% true if saving result as Zarr files
saveZarr = false;
% true if saving result as Uint16
Save16bit = true;

% use slurm cluster if true, otherwise use the local machine (master job)
parseCluster = false;
% use master job for task computing or not. 
masterCompute = true;
% configuration file for job submission
configFile = '';
% if true, use Matlab runtime (for the situation without matlab license)
mccMode = false;

XR_deskew_rotate_data_wrapper(dataPath_exps, xyPixelSize=xyPixelSize, dz=dz, ...
    Reverse=Reverse, largeFile=largeFile, zarrFile=zarrFile, saveZarr=saveZarr, ...
    Save16bit=Save16bit, parseCluster=parseCluster, masterCompute=masterCompute, ...
    configFile=configFile, mccMode=mccMode);


%% visualize raw, RL decon, and OMW decon

fsn = 'Scan_Iter_0000_0000_CamB_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_003y_000z_0000t';
fns = {[dataPath, 'DSR/', fsn, '.tif'], ...
       [dataPath,'matlab_decon_conventional/DSR/', fsn, '.tif'], ...
       [dataPath,'matlab_decon_omw/DSR/', fsn, '.tif'], ...       
      };

im_raw = readtiff(fns{1});
im_rl_decon = readtiff(fns{2});
im_omw_decon = readtiff(fns{3});

% crop a small region for visualization
bbox = [1040, 250, 133, 1240, 450, 133];

im_raw_crop = crop3d_mex(im_raw, bbox);
im_rl_decon_crop = crop3d_mex(im_rl_decon, bbox);
im_omw_decon_crop = crop3d_mex(im_omw_decon, bbox);
clear im_raw im_rl_decon im_omw_decon;

figure, 
set(gcf, 'color', 'w', 'position', [1, 1, 1200, 400])
t = tiledlayout(1, 3, "TileSpacing", 'compact', 'Padding', 'compact');
nexttile
imagesc(im_raw_crop);
clim([100, prctile(im_raw_crop(:), 99)]);
colormap('gray');
axis equal
axis off
title('Raw')

nexttile
imagesc(im_rl_decon_crop);
clim([prctile(im_rl_decon_crop(:), 40), prctile(im_rl_decon_crop(:), 99)]);
colormap('gray');
axis equal
axis off
title('RL 30 iters')

nexttile
imagesc(im_omw_decon_crop);
clim([prctile(im_omw_decon_crop(:), 40), prctile(im_omw_decon_crop(:), 99)]);
colormap('gray');
axis equal
axis off
title('OMW 2 iters')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compare conventional, BW, and OMW methods for the same image (skewd space)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fsn = 'Scan_Iter_0000_0000_CamB_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_003y_000z_0000t';
fn = [dataPath, fsn, '.tif'];
psfFn = [dataPath, 'PSF/488_2_c.tif'];

im = readtiff(fn);
psf = readtiff(psfFn);

% remove background of the image
bg = 100;
im = max(im - bg, 0);

% clean psf (remove background)

% data scan step size
dz_data = 0.3;
% psf scan step size
dz_psf = 0.3;
% median factor for thresholding, not used in the 'masked' method
medFactor = 1.5;
% psf generation method
PSFGenMethod = 'masked';
psf = psf_gen_new(psf, dz_psf, dz_data, medFactor, PSFGenMethod);
psf = psf ./ sum(psf, 'all');


%% conventional RL decon
% 30 iterations

nIter = 30;
tic
rl_decon = decon_lucy_function(im, psf, nIter);
toc


%% generate wb backprojector and decon
%
% the WB back projector is from the code in Guo et al. and included
% third_parties/WBDeconvolution

% back projector type
bp_type = 'wiener-butterworth';
% wiener alpha parameter, using suggested values from Guo et al. 2020. 
alpha = 0.0001;
% Butterworth filter parameter, using suggested values from Guo et al. 2020. 
beta = 0.001;
% order of the Butterworth filter, using suggested values from Guo et al. 2020. 
n = 10;
% resolution definition type, use FWHM for the demo
resFlag = 1;
% user defined resolution, we use the ones from PSF
iRes = [];
% visualization flag
verboseFlag = 1;
[psf_wb_b, ~] = BackProjector(psf, bp_type, alpha, beta, n, resFlag, iRes, verboseFlag);

% decon
nIter = 2;
tic
wb_decon = decon_lucy_omw_function(im, psf, psf_wb_b, nIter);
toc

%% omw decon

% OTF thresholding parameter
OTFCumThresh = 0.9;
% true if the PSF is in skew space
skewed = true;
XR_visualize_OTF_mask_segmentation(psfFn, OTFCumThresh, skewed);

% generate omw backprojector
% wiener filter parameter
alpha = 0.005;
% true if the PSF is in skew space
skewed = true;
% OTF cumulative percentile threshold for OTF mask segmentation
OTFCumThresh = 0.9;
% hann window range applied to the distance transform, 0.0 means the center and 1.0 means border of OTF mask
hanWinBounds = [0.8, 1.0];
[b_omw, OTF_bp_omw, abs_OTF_c, OTF_mask] = omw_backprojector_generation(psf, alpha, skewed, ...
    'OTFCumThresh', OTFCumThresh, 'hanWinBounds', hanWinBounds);

% decon
nIter = 2;
tic
omw_decon = decon_lucy_omw_function(im, psf, b_omw, nIter);
toc


%% deskew/rotate raw and deconvolved data and also crop a small region for visualization

% skew angle 
Angle = 32.45;
% scan step size
dz = 0.3;
% xy pixel size
xyPixelSize = 0.108;
% scan direction
Reverse = true;
% true if in objective scan mode, false for x scan mode
ObjectiveScan = false;
% resampling factor
resample = [];
% interpolation type
Interp = 'linear';

% bbox for crop
bbox = [1040, 250, 133, 1240, 450, 133];

raw_dsr = deskewRotateFrame3D(single(im), Angle, dz, xyPixelSize, ...
    'reverse', Reverse, 'Crop', true, 'ObjectiveScan', ObjectiveScan, ...
    'resample', resample, 'Interp', Interp);
raw_dsr_crop = crop3d_mex(uint16(raw_dsr), bbox);
clear raw_dsr;

rl_dsr = deskewRotateFrame3D(single(rl_decon), Angle, dz, xyPixelSize, ...
    'reverse', Reverse, 'Crop', true, 'ObjectiveScan', ObjectiveScan, ...
    'resample', resample, 'Interp', Interp);
rl_dsr_crop = crop3d_mex(uint16(rl_dsr), bbox);
clear rl_dsr;

wb_dsr = deskewRotateFrame3D(single(wb_decon), Angle, dz, xyPixelSize, ...
    'reverse', Reverse, 'Crop', true, 'ObjectiveScan', ObjectiveScan, ...
    'resample', resample, 'Interp', Interp);
wb_dsr_crop = crop3d_mex(uint16(wb_dsr), bbox);
clear wb_dsr;

omw_dsr = deskewRotateFrame3D(single(omw_decon), Angle, dz, xyPixelSize, ...
    'reverse', Reverse, 'Crop', true, 'ObjectiveScan', ObjectiveScan, ...
    'resample', resample, 'Interp', Interp);
omw_dsr_crop = crop3d_mex(uint16(omw_dsr), bbox);
clear omw_dsr;


%% visualize the cropped region


figure, 
set(gcf, 'color', 'w', 'position', [1, 1, 1600, 400])
t = tiledlayout(1, 4, "TileSpacing", 'compact', 'Padding', 'compact');

nexttile
imagesc(raw_dsr_crop);
clim([prctile(raw_dsr_crop(:), 1), prctile(raw_dsr_crop(:), 99)]);
colormap('gray');
axis equal
axis off
title('Raw')

nexttile
imagesc(rl_dsr_crop);
clim([prctile(rl_dsr_crop(:), 40), prctile(rl_dsr_crop(:), 99)]);
colormap('gray');
axis equal
axis off
title('RL 30 iters')

nexttile
imagesc(wb_dsr_crop);
clim([prctile(wb_dsr_crop(:), 40), prctile(wb_dsr_crop(:), 99)]);
colormap('gray');
axis equal
axis off
title('WB 2 iters')

nexttile
imagesc(omw_dsr_crop);
clim([prctile(omw_dsr_crop(:), 40), prctile(omw_dsr_crop(:), 99)]);
colormap('gray');
axis equal
axis off
title('OMW 2 iters')



