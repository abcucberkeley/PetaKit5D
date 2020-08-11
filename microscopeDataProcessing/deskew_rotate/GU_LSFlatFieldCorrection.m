function [Rawdata] = GU_LSFlatFieldCorrection(Rawdata,LSImage,background, varargin)

% Gokul Upadhyayula, 2016
% xruan (08/05/2020): add option for constant background after correction.

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('Rawdata'); %
ip.addRequired('LSImage'); %
ip.addRequired('background'); %;
ip.addParamValue('LowerLimit', 0.4, @isnumeric); % this value is the lowest
ip.addParameter('constOffset', [], @(x) isempty(x) || isnumeric(x)); % If it is set, use constant background, instead of background from the camera.
ip.parse(Rawdata, LSImage,background, varargin{:});
LowerLimit = ip.Results.LowerLimit;

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
LSImage = double(LSImage) - double(background);
LSImage = double(LSImage)/double(max(LSImage(:)));
Mask = repmat(LSImage,1,1,size(Rawdata,3));
Mask(Mask<LowerLimit) = LowerLimit;
background = repmat(background,1,1,size(Rawdata,3));

Temp = double(Rawdata);
Temp = Temp-double(background);
Temp(Temp<0) = 0;

% correct for odd numbered pixels
if size(Temp,1)<size(Mask,1)
    Mask(end,:,:) = [];
end
if size(Temp,2)<size(Mask,2)
    Mask(:,end,:) = [];
end

% LS Flat-field correction
if ~isempty(ip.Results.constOffset)
    Temp = Temp./double(Mask) + ip.Results.constOffset;
else
    Temp = Temp./double(Mask) + double(background);
end
Rawdata = single(Temp);
