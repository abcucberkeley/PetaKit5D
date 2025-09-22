function [relative_shift_yxz_mat] = XR_chromatic_shift_estimation(psfFullpaths, varargin)
% estimate chromatic shift with given PSFs


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('psfFullpaths', @(x) ischar(x) || iscell(x));
ip.addParameter('resultFullpath', '', @(x) ischar(x) || iscell(x));
ip.addParameter('maxOffset', [30, 30, 30], @isnumeric);
ip.addParameter('cropLength', [0, 0, 0], @isnumeric);
ip.addParameter('zarrFile', false, @islogical);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);

ip.parse(psfFullpaths, varargin{:});

pr = ip.Results;
resultFullpath = pr.resultFullpath;
maxOffset = pr.maxOffset;
cropLength = pr.cropLength;
zarrFile = pr.zarrFile;
uuid = pr.uuid;


nF = numel(psfFullpaths) ;
if nF < 2
    error('At least two PSF files are needed!');
end

% check relative shif aganist the first channel
if zarrFile
    im_s = readzarr(psfFullpaths{1});
else
    im_s = readtiff(psfFullpaths{1});
end

relative_shift_yxz_mat = zeros(nF, 3);
max_xcorr_mat = zeros(nF, 1);

for i = 2 : nF
    if zarrFile
        im_t = readzarr(psfFullpaths{i});
    else
        im_t = readtiff(psfFullpaths{i});
    end
    [relative_shift_yxz_i, max_xcorr] = chromatic_shift_estimation_for_two_channels(im_s, im_t, maxOffset, cropLength);
    relative_shift_yxz_mat(i, :) = relative_shift_yxz_i;
    max_xcorr_mat(i) = max_xcorr;
end

% save result file
if ~isempty(resultFullpath)
    if isempty(uuid)
        uuid = get_uuid();
    end
    [resultPath, fsn, ext] = fileparts(resultFullpath);
    resultTmpPath = sprintf('%s/%s_%s%s', resultPath, fsn, uuid, ext);
    save('-v7.3', resultTmpPath, 'psfFullpaths', 'maxOffset', 'cropLength', 'relative_shift_yxz_mat', 'max_xcorr_mat');
    movefile(resultTmpPath, resultFullpath);
end

end


function [relative_shift_yxz, max_xcorr] = chromatic_shift_estimation_for_two_channels(im_s, im_t, maxOffset, cropLength)
% use cross correlation registration
% axis order for all parameters are in yxz

im_t1 = im_t(1 + cropLength(1) : end - cropLength(1), 1 + cropLength(2) : end - cropLength(2), 1 + cropLength(3) : end - cropLength(3));

maxShifts_lb = cropLength - maxOffset;
maxShifts_ub = cropLength + maxOffset;
maxShifts = [maxShifts_lb; maxShifts_ub];
maxShifts = maxShifts(:, [2, 1, 3]);

[offset_yxz, max_xcorr, C] = normxcorr3_max_shift(single(im_t1), single(im_s), maxShifts);
sp2 = -cropLength;
relative_shift_yxz = offset_yxz + sp2; 

end