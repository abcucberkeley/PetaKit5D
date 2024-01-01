function writetiff(img, filepath, options)
% wrapper for saveastiff.m to support big tiff format, and keep
% compatibility of previous script for the decrapted writetiff. 
% 
% Author: Xiongtao Ruan (07/30/2020)
% xruan (08/07/2020): add option for group write
% xruan (03/24/21): fix issue for compession method not in lower case
% xruan (07/21/2022): add parallelWriteTiff as default method


arguments
    img {mustBeNumeric}
    filepath char
    options.Compression (1, :) char {mustBeMember(options.Compression, {'none', 'lzw'})} = 'lzw'
    options.Mode (1, :) char {mustBeMember(options.Mode, {'parallel', 'libtiff', 'imwrite'})} = 'parallel'
    options.groupWrite (1, 1) logical = true
end

Compression = options.Compression;
Mode = options.Mode;
groupWrite = options.groupWrite;

options = struct();
switch lower(Compression)
    case 'none'
        options.compress = 'no';
    case 'lzw'
        options.compress = 'lzw';
end

switch lower(Mode)
    case 'parallel'
        try
            pstr = fileparts(filepath);
            if ~exist(pstr, 'dir')
                mkdir_recursive(pstr);
            end
            parallelWriteTiff(filepath, img, 'w');
        catch ME
            disp(ME)
            disp('Use the alternative tiff writer (saveastiff)...');
            options.message = false;
            options.overwrite = true;
            saveastiff(img, filepath, options);            
        end
    case 'libtiff'
        options.message = false;
        options.overwrite = true;
        saveastiff(img, filepath, options);
end

if groupWrite && ~ispc
    try
        fileattrib(filepath, '+w', 'g');
    catch
        warning('Unable to change file attribe for group write!');
    end
end

end

