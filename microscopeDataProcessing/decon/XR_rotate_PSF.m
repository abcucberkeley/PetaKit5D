function [] = XR_rotate_PSF(PSFfile, varargin)
% function to rotate PSF, especally for the rotated data with sample scan.
% 
% Author: Xiongtao Ruan (07/12/2020)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('PSFfile'); 
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('ObjectiveScan', false, @islogical);
ip.addParameter('xyPixelSize', 0.108, @isscalar);
ip.addParameter('dz', 0.1, @isscalar);
ip.addParameter('SkewAngle', 32.45, @isscalar);
ip.addParameter('HighPrctile', 99, @isscalar);
ip.addParameter('Reverse', false, @islogical);
ip.addParameter('Save16bit', true , @islogical); % saves deskewed data as 16 bit -- not for quantification
ip.addParameter('uuid', '', @isstr);

ip.parse(PSFfile, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
xyPixelSize = pr.xyPixelSize;
dz = pr.dz;
Reverse = pr.Reverse;
SkewAngle = pr.SkewAngle;
ObjectiveScan = pr.ObjectiveScan;
Save16bit = pr.Save16bit;

uuid = ip.Results.uuid;
% uuid for the job
if isempty(uuid)
    uuid = get_uuid();
end

[PSFdir, fsname] = fileparts(PSFfile);
rotPSFdir = [PSFdir, '/', 'Rotated/'];
if ~exist(rotPSFdir, 'dir')
    mkdir(rotPSFdir);
end
rtFullname = sprintf('%s%s.tif', rotPSFdir, fsname);
if exist(rtFullname, 'file') && ~Overwrite
    return;
end

% decide zAniso
if ObjectiveScan
    zAniso = dz / xyPixelSize;
else
    theta = SkewAngle * pi / 180;
    dz0 = sin(theta) * dz; % ~0.25 for dz0 = 0.45
    zAniso = dz0 / xyPixelSize;
end

im = double(readtiff(PSFfile));
im_rt = rotateFrame3D(im, SkewAngle, zAniso, Reverse, 'Crop', true, 'ObjectiveScan', ObjectiveScan);

% pad 0 with median 
med = median(im_rt(im_rt > 0 & im_rt < prctile(im_rt(im_rt > 0), pr.HighPrctile)));
im_rt(im_rt == 0) = med;

% save the rotated PSF files to the subfolder of PSF folder and also save
% the parameters. 

rtTempName = sprintf('%s%s_%s.tif', rotPSFdir, fsname, uuid);
if Save16bit
    writetiff(uint16(im_rt), rtTempName);
else
    writetiff(single(im_rt), rtTempName);
end
movefile(rtTempName, rtFullname);

% save parameters 
save('-v7.3', [rotPSFdir, '/parameters.mat'], 'pr');
writetable(struct2table(pr, 'AsArray', true), [rotPSFdir, '/parameters.txt'])

end
