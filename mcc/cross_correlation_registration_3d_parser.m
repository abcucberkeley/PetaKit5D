function [] = cross_correlation_registration_3d_parser(imgFullpath_1, imgFullpath_2, xcorrFullpath, cuboid_1, cuboid_2, cuboid_overlap_12, px, xyz_factors, varargin)


%#function cross_correlation_registration_2d

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('imgFullpath_1', @ischar);
ip.addRequired('imgFullpath_2', @ischar);
ip.addRequired('xcorrFullpath', @ischar);
ip.addRequired('cuboid_1', @(x) isnumeric(x) || ischar(x));
ip.addRequired('cuboid_2', @(x) isnumeric(x) || ischar(x));
ip.addRequired('cuboid_overlap_12', @(x) isnumeric(x) || ischar(x));
ip.addRequired('px', @(x) isnumeric(x) || ischar(x));
ip.addRequired('xyz_factors', @(x) isnumeric(x) || ischar(x));
ip.addParameter('downSample', [1, 1, 1], @(x) isnumeric(x) || ischar(x));
ip.addParameter('MaxOffset', [300, 300, 50], @(x) isnumeric(x) || ischar(x));
ip.addParameter('dimNumThrsh', 10000, @(x) isnumeric(x) || ischar(x));

ip.parse(imgFullpath_1, imgFullpath_2, xcorrFullpath, cuboid_1, cuboid_2, cuboid_overlap_12, px, xyz_factors, varargin{:});

pr = ip.Results;
downSample = pr.downSample;
MaxOffset = pr.MaxOffset;
dimNumThrsh = pr.dimNumThrsh;

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
if ischar(downSample)
    downSample = str2num(downSample);
end
if ischar(MaxOffset)
    MaxOffset = str2num(MaxOffset);
end
if ischar(dimNumThrsh)
    dimNumThrsh = str2double(dimNumThrsh);
end

cross_correlation_registration_3d(imgFullpath_1, imgFullpath_2, xcorrFullpath, ...
    cuboid_1, cuboid_2, cuboid_overlap_12, px, xyz_factors, downSample=downSample, ...
    MaxOffset=MaxOffset, dimNumThrsh=dimNumThrsh);


end
