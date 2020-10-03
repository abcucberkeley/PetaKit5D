function [imSize] = getImageSize(filePath)
% Get image size without load the image
% 
% use Tiff for faster performance

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
% toc


end

