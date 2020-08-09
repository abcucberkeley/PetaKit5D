function [ResidualSigmaThresh] = XR_get_camera_standard_deviation_threshold(data, varargin)
% The function is to set background threshold for standard deviation based
% on camera
% usually base level std for SMOS is 4, and it is ~54 for EMCCD. 
% The idea is to obtain the camera information and set the threshold. If
% the camera information is not avaiable, estimate the threshold from
% image. 
% 
% Author: Xiongtao Ruan 11/08/2019
% 


ip = inputParser;
ip.addRequired('data', @isstruct);
ip.addParameter('SCMOS_DC_Thresh', 3, @isnumeric);
ip.addParameter('EMCCD_DC_Thres', 65, @isnumeric);
ip.parse(data, varargin{:});

camera_class = get_camera_type_from_setting_file(data);

switch camera_class
    case 'SCMOS'
        ResidualSigmaThresh = ip.Results.SCMOS_DC_Thresh;
    case 'EMCCD'
        ResidualSigmaThresh = ip.Results.EMCCD_DC_Thres;
    case ''
       ResidualSigmaThresh = nan;
end


end
