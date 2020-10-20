function [imSize] = getImageSize(filePath)
% Get image size without load the image
% 
% use Tiff for faster performance
% 
% Author: Xiongtao Ruan (01/2020)
% xruan (10/15/2020): add support for zarr format

warning ('off','all');

% tic
% fileInfo = imfinfo(filePath);
% 
% Width = fileInfo(1).Width;
% Height = fileInfo(1).Height;
% Zstack = numel(fileInfo);
% 
% imSize = [Height, Width, Zstack];
% toc

% tic

if strcmp(filePath(end - 2 : end), 'tif') || strcmp(filePath(end - 3 : end), 'tiff')
    t = Tiff(filePath, 'r');
    Height = getTag(t, 'ImageLength');
    Width = getTag(t, 'ImageWidth');
    Zstack = 1;
    while ~t.lastDirectory()
        t.nextDirectory() ;
        Zstack = Zstack + 1;
    end
    close(t);
    imSize = [Height, Width, Zstack];
elseif strcmp(filePath(end - 3 : end), 'zarr')
    bim = blockedImage(filePath, 'Adapter', ZarrAdapter);
    imSize = bim.Size;
end    
    
    
% toc

% tic
% t = matlab.io.internal.BigImageTiffReader(filePath);
% imSize = [t.ImageHeight, t.ImageWidth, t.NumImages];
% toc



end

