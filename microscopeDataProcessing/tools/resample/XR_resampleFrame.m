function [] = XR_resampleFrame(fn, fnout, rsfactor, varargin)
% resample image for a single frame
%
% 
% Author: Xiongtao Ruan (01/15/2021)
% 
% xruan (05/11/2022): add support for 2d image

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('fn', @ischar);
ip.addRequired('fnout', @ischar);
ip.addRequired('rsfactor', @isnumeric);
ip.addParameter('Interp', 'linear', @ischar);
ip.addParameter('Save16bit', true ,@islogical); % saves 16bit, else single
ip.addParameter('zarrFile', false, @islogical);
ip.addParameter('saveZarr', false, @islogical); % use zarr file as output
ip.parse(fn, fnout, rsfactor, varargin{:});

ip.parse(fn, fnout, rsfactor, varargin{:});
pr = ip.Results;

Interp = pr.Interp;
Save16bit = pr.Save16bit;
zarrFile = pr.zarrFile;
saveZarr = pr.saveZarr;

if exist(fnout, 'file')
    fprintf('The result file %s already exists, skip it!', fnout);
    return;
end

if zarrFile
    im = readzarr(fn);
else
    try 
        im = parallelReadTiff(fn);
    catch ME
        disp(ME);
        im = readtiff(fn);
    end
end
im = single(im);

rs = rsfactor(:)';
if ismatrix(im)
    rs = [ones(1, 3 - numel(rs)) * rs(1), rs(2:end)];    
else
    % complete rs to 3d in case it is not
    rs = [ones(1, 4 - numel(rs)) * rs(1), rs(2:end)];    
end

outSize = round(size(im) ./ rs);

if ismatrix(im)
    im = imresize(im, outSize, 'Method', Interp);    
else
    im = imresize3(im, outSize, 'Method', Interp);
end

if Save16bit
    im = uint16(im);
end

[pstr, fsname, ext] = fileparts(fnout);
uuid = get_uuid();
fnout_tmp = [pstr, '/', fsname, '_', uuid, ext];
if saveZarr
    writezarr(im, fnout_tmp);
else
    writetiff(im, fnout_tmp);
end
movefile(fnout_tmp, fnout);

end

