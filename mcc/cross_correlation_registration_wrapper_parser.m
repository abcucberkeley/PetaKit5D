function [] = cross_correlation_registration_wrapper_parser(imgFullpath_1, imgFullpath_2, xcorrFullpath, pair_indices, cuboid_1, cuboid_2, cuboid_overlap_12, px, xyz_factors, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('imgFullpath_1', @ischar);
ip.addRequired('imgFullpath_2', @(x) iscell(x) || ischar(x));
ip.addRequired('xcorrFullpath', @ischar);
ip.addRequired('pair_indices', @(x) isnumeric(x) || ischar(x));
ip.addRequired('cuboid_1', @(x) isnumeric(x) || ischar(x));
ip.addRequired('cuboid_2', @(x) isnumeric(x) || ischar(x));
ip.addRequired('cuboid_overlap_12', @(x) isnumeric(x) || ischar(x));
ip.addRequired('px', @(x) isnumeric(x) || ischar(x));
ip.addRequired('xyz_factors', @(x) isnumeric(x) || ischar(x));
ip.addParameter('Stitch2D', false, @(x) islogical(x) || ischar(x));
ip.addParameter('downSample', [1, 1, 1], @(x) isnumeric(x) || ischar(x));
ip.addParameter('MaxOffset', [300, 300, 50], @(x) isnumeric(x) || ischar(x));
ip.addParameter('largeZarr', false, @(x) islogical(x) || ischar(x));
ip.addParameter('mipDirStr', '', @ischar);
ip.addParameter('poolSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('dimNumThrsh', 10000, @(x) isnumeric(x) || ischar(x));

ip.parse(imgFullpath_1, imgFullpath_2, xcorrFullpath, pair_indices, cuboid_1, cuboid_2, cuboid_overlap_12, px, xyz_factors, varargin{:});

pr = ip.Results;
Stitch2D = pr.Stitch2D;
downSample = pr.downSample;
MaxOffset = pr.MaxOffset;
largeZarr = pr.largeZarr;
mipDirStr = pr.mipDirStr;
poolSize = pr.poolSize;
dimNumThrsh = pr.dimNumThrsh;

if ischar(imgFullpath_2)
    imgFullpath_2 = eval(imgFullpath_2);
end
if ischar(pair_indices)
    pair_indices = str2num(pair_indices);
end
if ischar(cuboid_1)
    cuboid_1 = str2num(cuboid_1);
end
if ischar(cuboid_2)
    cuboid_2 = str2num(cuboid_2);
end
if ischar(cuboid_overlap_12)
    cuboid_overlap_12 = str2num(cuboid_overlap_12);
end
if ischar(px)
    px = str2double(px);
end
if ischar(xyz_factors)
    xyz_factors = str2num(xyz_factors);
end
if ischar(Stitch2D)
    Stitch2D = strcmp(Stitch2D, 'true');
end
if ischar(downSample)
    downSample = str2num(downSample);
end
if ischar(MaxOffset)
    MaxOffset = str2num(MaxOffset);
end
if ischar(largeZarr)
    largeZarr = strcmp(largeZarr, 'true');
end
if ischar(poolSize)
    poolSize = str2num(poolSize);
end
if ischar(dimNumThrsh)
    dimNumThrsh = str2double(dimNumThrsh);
end

cross_correlation_registration_wrapper(imgFullpath_1, imgFullpath_2, xcorrFullpath, ...
    pair_indices, cuboid_1, cuboid_2, cuboid_overlap_12, px, xyz_factors, Stitch2D=Stitch2D, ...
    downSample=downSample, MaxOffset=MaxOffset, largeZarr=largeZarr, mipDirStr=mipDirStr, ...
    poolSize=poolSize, dimNumThrsh=dimNumThrsh);

end
