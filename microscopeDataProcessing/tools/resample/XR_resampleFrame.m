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
ip.addParameter('bbox', [], @isnumeric); % bbox for input
ip.addParameter('Interp', 'linear', @ischar);
ip.addParameter('Save16bit', true ,@islogical); % saves 16bit, else single
ip.addParameter('zarrFile', false, @islogical);
ip.addParameter('saveZarr', false, @islogical); % use zarr file as output
ip.addParameter('blockSize', [256, 256, 256], @isnumeric); % blcoksize
ip.addParameter('uuid', '', @ischar);

ip.parse(fn, fnout, rsfactor, varargin{:});

pr = ip.Results;
bbox = pr.bbox;
Interp = pr.Interp;
Save16bit = pr.Save16bit;
zarrFile = pr.zarrFile;
saveZarr = pr.saveZarr;
blockSize = pr.blockSize;
uuid = pr.uuid;

if isempty(uuid)
    uuid = get_uuid();
end

if exist(fnout, 'file')
    fprintf('The result file %s already exists, skip it!', fnout);
    return;
end

if zarrFile
    im = readzarr(fn, bbox=bbox);
else
    if isempty(bbox)
        im = readtiff(fn);
    else
        im = readtiff(fn, range=[bbox(3), bbox(6)]);
        bbox_1 = [bbox(1), bbox(2), 1, bbox(4), bbox(5), bbox(6) - bbox(3) + 1];
        try
            im = crop3d_mex(im, bbox_1);
        catch ME
            disp(ME)
            fprintf('Use alternative crop method...\n');
            im = im(bbox_1(1) : bbox_1(4), bbox_1(2) : bbox_1(5), :);
        end
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
    if strcmpi(Interp, 'linear')
        Interp = 'bilinear';
    end
    im = imresize(im, outSize, 'Method', Interp);    
else
    im = imresize3(im, outSize, 'Method', Interp);
end

if Save16bit
    im = uint16(im);
end

[pstr, fsname, ext] = fileparts(fnout);
fnout_tmp = [pstr, '/', fsname, '_', uuid, ext];
if saveZarr
    writezarr(im, fnout_tmp, blockSize=blockSize);
else
    writetiff(im, fnout_tmp);
end
movefile(fnout_tmp, fnout);

end

