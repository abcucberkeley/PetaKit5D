% demo to analysis data from scratch with the integrated function. 

%% Step 1: load data
% add the software to the path
addpath(genpath('../'));

%  rt is the root directory of data, you may provide root directory. If it
%  is not provided, the load condition data function will prompt out to
%  choose the root directory interactively. 
rt = '/clusterfs/fiona/xruan/Images/LLSM3DTools_test_images/2015-04-7_p505_p6_sumCLTArfp_shiga/Ex1_2minAfter_560_50mW_20p_647_150mW_20p_40ms_z0p5';

% follow the below when asking for inputs (They may be different for your own experiments):
% Enter the number of channels: 2
% click folder name 'ch1' when it first prompts out. 
% click folder name 'ch2' when it prompts out again. 
% Enter the fluorescent marker for channel 1: gfp
% Enter the fluorescent marker for channel 2: rfp
data = XR_loadConditionData3D(rt);


%% Step 2 estimate psf sigmas if there are calibration files and they are not available (optional)
% The sigmas of psfs are estimated separately. The filename is provided as 
% input for the estimation. 
ch1_psf_filename = '/clusterfs/fiona/xruan/Images/LLSM3DTools_test_images/2015-04-7_p505_p6_sumCLTArfp_shiga/560totalPSF.tif';
[sigmaXY_ch1, sigmaZ_ch1] = GU_estimateSigma3D(ch1_psf_filename, []);
ch2_psf_filename = '/clusterfs/fiona/xruan/Images/LLSM3DTools_test_images/2015-04-7_p505_p6_sumCLTArfp_shiga/642totalPSF.tif';
[sigmaXY_ch2, sigmaZ_ch2] = GU_estimateSigma3D(ch2_psf_filename, []);

% Because in calibration (of this data), the calibration image is
% approximately isotropic, which is is not true for the data. We need to
% divide the sigma in z-axis by the Anisotropic factor in this axis. 
sigma_mat = [sigmaXY_ch1, sigmaZ_ch1 ./ data(1).zAniso; 
             sigmaXY_ch2, sigmaZ_ch2 ./ data(1).zAniso];


%% step 3: deskew, detection, tracking and tracking postprocessing
% The function is the main script that integrate deskew, detection
% tracking, tracking postprocessing, and tracking rotation

XR_runDetTrack3d(data, 'Sigma', sigma_mat, ..., 
    'DetectionMethod', 'lowSNR', 'aname', 'Analysis_update_3', ...
    'FitMixtures', false, ...
    'BackgroundCorrection', true, 'XZoffsetCorrection', true);
