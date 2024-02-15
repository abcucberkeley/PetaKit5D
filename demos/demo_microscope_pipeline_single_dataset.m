% demo to run microscope automatic pipeline for a single dataset with stitch

% note: this demo is designed for quick processing of raw microscopy data 
% during acquistion, providing fast feedback and inspection. For final
% processing after acquistion, follow the functions outlined in the stitching, 
% deconvolution, geometric transformation, and large-scale processing demos.


clear, clc;

fprintf('Microscopy pipeline for a single dataset...\n\n');

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


%% Step 2: set parameters for each step. 
% add the software to the path
cd(fileparts(which(mfilename)));
addpath(genpath('../'));

% data path
% defined above

% general parameters.
% pixel size in xy-plane
xyPixelSize = 0.108;
% pixel size in z-axis
dz = 0.3;
% Inverse direction of z axis. 
Reverse = true;
% Channel identifiers and orders. 
ChannelPatterns = {'CamB_ch0', 'CamA_ch0'};
% Save 16bit for deskew/rotate, stitch, deconvolution, and rotation after decon
Save16bit = [true, true, true, true];
% not overwrite existing results
Overwrite = false;
% Processing for existing dataset
Streaming = false;
% number of cores for the pipeline (only for ABC cluster). For other
% environment, this setting will be ignored. 
cpusPerTask = 24;

% define which step is included in the analysis
Deskew = true; % deskew
Rotate = true; % rotate after deskew
Stitch = true; % stitch
Decon = false; % deconvolution
RotateAfterDecon = false; % rotate after deconvolution


% flat field parameters (if LLFFCorrection is set to true) 
LLFFCorrection = true;
% flat field image paths in the order of CamB_ch0 and CamA_ch0
LSImagePaths = {[dataPath, '/FF/averaged/ff_CamB_ch0_CAM1_stack0000_488nm_0000000msec_0096911533msecAbs_000x_000y_000z_0017t.tif'],...
                [dataPath, '/FF/averaged/ff_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0096911533msecAbs_000x_000y_000z_0017t.tif']};

% camera background image paths in the order of CamB and CamA
BackgroundPaths = {[dataPath, '/FF/KorraFusions/AVG_DF400_CamB_10ms.tif'],...
                   [dataPath, '/FF/KorraFusions/AVG_DF400_CamA_10ms.tif']};

% stitch parameters
% image list csv file path
ImageListFullpaths = [dataPath, 'ImageList_from_encoder.csv'];
% axis order for coordinates
axisOrder = '-x,y,z';


%% Step 3: run the analysis with given parameters. 

XR_microscopeAutomaticProcessing(dataPath, 'xyPixelSize', xyPixelSize, 'dz', dz,  ...
    'Reverse', Reverse, 'ChannelPatterns', ChannelPatterns, 'Save16bit', Save16bit, ...
    'Overwrite', Overwrite, 'Streaming', Streaming, 'Deskew', Deskew, 'Rotate', Rotate, ...
    'Stitch', Stitch, 'Decon', Decon, 'RotateAfterDecon', RotateAfterDecon, ...
    'ImageListFullpaths', ImageListFullpaths, 'axisOrder', axisOrder, 'cpusPerTask', cpusPerTask);

