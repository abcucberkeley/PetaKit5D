function [dtype] = getImageDataType(filePath, zarrFile)
% Get image data type without loading the image
% 
% 
% Author: Xiongtao Ruan (04/2023)


% warning ('off','all');

if nargin < 2
    zarrFile = false;
end

if strcmp(filePath(end - 2 : end), 'tif') || strcmp(filePath(end - 3 : end), 'tiff')
    % directly call Big Image Tiff Reader using in MPageTiffAdapter.m
    tobj = matlab.io.internal.BigImageTiffReader(filePath);
    dtype = tobj.MLType;
elseif strcmp(filePath(end - 3 : end), 'zarr') || zarrFile
    try 
        bim = blockedImage(filePath, 'Adapter', CZarrAdapter);
    catch ME
        disp(ME);
        bim = blockedImage(filePath, 'Adapter', ZarrAdapter);
    end
    dtype = bim.ClassUnderlying;
else
    error('Unknown data type, currently only tiff and zarr files are supported!')
end


end

