function [Temp] = XR_LSFlatFieldCorrection(Rawdata, LSImage,background, varargin)

% Gokul Upadhyayula, 2016
% xruan (08/05/2020): add option for constant background after correction.
% xruan (09/24/2021): keep LSImage and background as 2d and use
% broadcasting for the computing to reduce required memory. Also change to
% use single format to save memory.
% xruan (07/15/2022): add support to not subtract background from LS Image
% (because in a lot of situations, the LS Image is estimated with no background).
% xruan (05/26/2023): rename it in case of conflict with old versions

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('Rawdata'); %
ip.addRequired('LSImage'); %
ip.addRequired('background'); %;
ip.addParameter('LowerLimit', 0.4, @isnumeric); % this value is the lowest
ip.addParameter('constOffset', [], @(x) isempty(x) || isnumeric(x)); % If it is set, use constant background, instead of background from the camera.
ip.addParameter('removeFFImBackground', true, @(x) islogical(x)); % true: subtract background; false: not subtract background
ip.addParameter('LSRescale', true, @(x) islogical(x)); % true: rescale LS by maximum; false: use flat field as it is. 
ip.addParameter('castDataType', true, @(x) islogical(x)); % true: cast the data type the same as input, otherwise output as single. 
ip.parse(Rawdata, LSImage,background, varargin{:});

pr = ip.Results;
LowerLimit = pr.LowerLimit;
constOffset = pr.constOffset;
removeFFImBackground = pr.removeFFImBackground;
LSRescale = pr.LSRescale;
castDataType = pr.castDataType;

% average z planes of LS image
if ndims(LSImage)==3
    LSImage = squeeze(mean(LSImage,3));
end

% crop LS data if necessary
ImSize = size(squeeze(Rawdata(:,:,1)));
MapSize = size(LSImage);
D =ceil((MapSize-ImSize)/2);
LSImage = LSImage(D(1)+1:D(1)+ImSize(1),D(2)+1:D(2)+ImSize(2));

% crop Bk data if necessary
MapSize = size(background);
D =ceil((MapSize-ImSize)/2);
background = background(D(1)+1:D(1)+ImSize(1),D(2)+1:D(2)+ImSize(2));

% Prepare LS Flat-field correction mask
if removeFFImBackground
    LSImage = single(LSImage) - single(background);
end
if LSRescale
    LSImage = single(LSImage)/single(max(LSImage(:)));
end
% Mask(Mask<LowerLimit) = LowerLimit;
LSImage = max(LSImage, LowerLimit);
% background = repmat(background,1,1,size(Rawdata,3));

try 
    LSImage = single(LSImage);
    background = single(background);
    if ~isempty(constOffset)
        Temp = flat_field_correction_3d_mex(Rawdata, LSImage, background, castDataType, constOffset);
    else
        Temp = flat_field_correction_3d_mex(Rawdata, LSImage, background, castDataType);        
    end
catch ME
    disp(ME)
    % Temp = double(Rawdata);
    Temp = single(Rawdata);
    dtype = class(Rawdata);
    clear Rawdata;
    Temp = Temp-single(background);
    % Temp(Temp<0) = 0;
    % Temp = Temp .* (Temp >= 0);
    Temp = max(Temp, 0);
    
    % LS Flat-field correction
    if ~isempty(constOffset)
        Temp = Temp./single(LSImage) + constOffset;
    else
        Temp = Temp./single(LSImage) + single(background);
    end
    if castDataType 
        Temp = cast(Temp, dtype);
    end
end

end
