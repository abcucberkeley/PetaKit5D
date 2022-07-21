function writetiff(img, filepath, varargin)
% wrapper for saveastiff.m to support big tiff format, and keep
% compatibility of previous script for the decrapted writetiff. 
% 
% Author: Xiongtao Ruan (07/30/2020)
% xruan (08/07/2020): add option for group write
% xruan (03/24/21): fix issue for compession method not in lower case
% xruan (07/21/2022): add parallelWriteTiff as default method


ip = inputParser;
ip.CaseSensitive = false;
ip.addParameter('Compression', 'lzw', @(x) any(strcmpi(x, {'none', 'lzw'})));
ip.addParameter('Mode', 'parallel', @(x) any(strcmpi(x, {'parallel', 'libtiff', 'imwrite'})));
ip.addParameter('groupWrite', true, @islogical);
ip.parse(varargin{:});

pr = ip.Results;

options = struct();
switch lower(pr.Compression)
    case 'none'
        options.compress = 'no';
    case 'lzw'
        options.compress = 'lzw';
end

switch lower(pr.Mode)
    case 'parallel'
        try
            parallelWriteTiff(filepath, img, 'w');
        catch ME
            disp(ME)
            options.message = false;
            options.overwrite = true;
            saveastiff(img, filepath, options);            
        end
    case 'libtiff'
        options.message = false;
        options.overwrite = true;
        saveastiff(img, filepath, options);
end

if pr.groupWrite
    try
        fileattrib(filepath, '+w', 'g');
    catch
        warning('Unable to change file attribe for group write!');
    end
end

end

