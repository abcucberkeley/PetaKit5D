function [imSize] = getImageSize(filePath)
% Get image size without load the image

warning ('off','all');

fileInfo = imfinfo(filePath);

Width = fileInfo(1).Width;
Height = fileInfo(1).Height;
Zstack = numel(fileInfo);

imSize = [Height, Width, Zstack];

end

