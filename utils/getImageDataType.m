function [dtype] = getImageDataType(filePath)
% Get image data type without loading the image
% 
% 
% Author: Xiongtao Ruan (04/2023)


% warning ('off','all');

if strcmp(filePath(end - 2 : end), 'tif') || strcmp(filePath(end - 3 : end), 'tiff')
    bim = blockedImage(filePath, 'Adapter', MPageTiffAdapter);
elseif strcmp(filePath(end - 3 : end), 'zarr')
    try 
        bim = blockedImage(filePath, 'Adapter', CZarrAdapter);
    catch ME
        disp(ME);
        bim = blockedImage(filePath, 'Adapter', ZarrAdapter);
    end
end

dtype = bim.ClassUnderlying;


end

