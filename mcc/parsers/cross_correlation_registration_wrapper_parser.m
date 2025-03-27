function [] = cross_correlation_registration_wrapper_parser(imgFullpath_1, imgFullpath_2, xcorrFullpath, pair_indices, cuboid_1, cuboid_2, cuboid_overlap_12, xyz_voxelsizes, data_order_mat, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('imgFullpath_1', @ischar);
ip.addRequired('imgFullpath_2', @(x) iscell(x) || ischar(x));
ip.addRequired('xcorrFullpath', @ischar);
ip.addRequired('pair_indices', @(x) isnumeric(x) || ischar(x));
ip.addRequired('cuboid_1', @(x) isnumeric(x) || ischar(x));
ip.addRequired('cuboid_2', @(x) isnumeric(x) || ischar(x));
ip.addRequired('cuboid_overlap_12', @(x) isnumeric(x) || ischar(x));
ip.addRequired('xyz_voxelsizes', @(x) isnumeric(x) || ischar(x));
ip.addRequired('data_order_mat', @(x) isnumeric(x) || ischar(x));
ip.addParameter('Stitch2D', false, @(x) islogical(x) || ischar(x));
ip.addParameter('downSample', [1, 1, 1], @(x) isnumeric(x) || ischar(x));
ip.addParameter('MaxOffset', [300, 300, 50], @(x) isnumeric(x) || ischar(x));
ip.addParameter('largeFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('mipDirStr', '', @ischar);
ip.addParameter('poolSize', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('dimNumThrsh', 10000, @(x) isnumeric(x) || ischar(x));
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('configFile', '', @ischar);

ip.parse(imgFullpath_1, imgFullpath_2, xcorrFullpath, pair_indices, cuboid_1, cuboid_2, cuboid_overlap_12, xyz_voxelsizes, data_order_mat, varargin{:});

pr = ip.Results;
Stitch2D = pr.Stitch2D;
downSample = pr.downSample;
MaxOffset = pr.MaxOffset;
largeFile = pr.largeFile;
mipDirStr = pr.mipDirStr;
poolSize = pr.poolSize;
dimNumThrsh = pr.dimNumThrsh;
parseCluster = pr.parseCluster;
mccMode = pr.mccMode;
configFile = pr.configFile;

if ischar(imgFullpath_2) && ~isempty(imgFullpath_2) && strcmp(imgFullpath_2(1), '{')
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
if ischar(xyz_voxelsizes)
    xyz_voxelsizes = str2num(xyz_voxelsizes);
end
if ischar(data_order_mat)
    data_order_mat = str2num(data_order_mat);
end
if ischar(Stitch2D)
    Stitch2D = str2num(Stitch2D);
end
if ischar(downSample)
    downSample = str2num(downSample);
end
if ischar(MaxOffset)
    MaxOffset = str2num(MaxOffset);
end
if ischar(largeFile)
    largeFile = str2num(largeFile);
end
if ischar(poolSize)
    poolSize = str2num(poolSize);
end
if ischar(dimNumThrsh)
    dimNumThrsh = str2num(dimNumThrsh);
end
if ischar(parseCluster)
    parseCluster = str2num(parseCluster);
end
if ischar(mccMode)
    mccMode = str2num(mccMode);
end

cross_correlation_registration_wrapper(imgFullpath_1, imgFullpath_2, xcorrFullpath, ...
    pair_indices, cuboid_1, cuboid_2, cuboid_overlap_12, xyz_voxelsizes, data_order_mat, ...
    Stitch2D=Stitch2D, downSample=downSample, MaxOffset=MaxOffset, largeFile=largeFile, ...
    mipDirStr=mipDirStr, poolSize=poolSize, dimNumThrsh=dimNumThrsh, parseCluster=parseCluster, ...
    mccMode=mccMode, configFile=configFile);

end

