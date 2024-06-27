function [] = XR_generate_image_list_wrapper(dataPaths, generationMethod, varargin)
% Wrapper to generate image list for a given generation method: encoder, sqlite, or tile_position 
% It saves as a csv file with name ImageList_from_*.csv in the dataPaths. The format is consistent with old csv files. 
% 
% 
% Required inputs:
%           dataPaths : char or cell array. Directory paths for the datasets. Either a string for a single dataset or a cell array of paths for several datasets with same experimental settings.
%    generationMethod : 'encoder'|'sqlite'|'tile_position'. Image list generation method. 'encoder': from encoder positions; 'sqlite': from sqlite database; 'tile_position': estimated from user-provided overlap size between neighboring tiles.
%
% Parameters (as 'specifier'-value pairs):
%     channelPatterns : a cell array (default: {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}).  Channel identifiers for included channels. 
%        tilePatterns : a 1x5 cell array (default: {'0000t', 'ch0', '000x', '000y', '000z'}). Patterns for time, channel, x, y and z to localize tiles. It should be the combination of word and numbers in the form of [a-zA-Z]*[0-9]* or [0-9]*[a-zA-Z]*.
%                  DS : true|false (default: false). Data is in deskewed space.
%                 DSR : true|false (default: false). Data is in deskew/rotated space (with stage coordinates).
%         xyPixelSize : a number (default: 0.108). Pixel size in um.
%                  dz : a number (default: 0.5). Scan interval in um.
%           skewAngle : a number (default: 32.45). Skew angle (in degree) of the stage.
%           axisOrder : char (default: 'xyz'). Axis order mapping for coordinates in image list. With combinations of -, x, y, and z. '-yxz' means negative -y map to x, x maps to y, and z maps to z.
%           dataOrder : char (default: 'y,x,z'). Axis order mapping for data. 'y,x,z' means the first, second and third axes are y, x, and z, respectively.
%       objectiveScan : true|false (default: false). Objective scan.
%              IOScan : true|false (default: false). Inverted objective scan. This is the scan with the stage coordinates (DSR space). 
%            zarrFile : true|false (default: false). Use Zarr file as input.
%         overlapSize : empty or 1x3 vector (default: []). Overlap size between tiles in pixel or um. If in pixels, axis order is yxz; if in um, axis order is xyz.
%     overlapSizeType : 'pixel'|'um' (default: 'pixel'). The unit for the overlap size.
%                uuid : empty or a uuid string (default: ''). uuid string as part of the temporate result paths.
%
%
% Author: Xiongtao Ruan (06/07/2024)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addRequired('generationMethod', @(x) ischar(x));
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @iscell);
ip.addParameter('tilePatterns', {'0000t', 'ch0', '000x', '000y', '000z'}, @iscell);
ip.addParameter('DS', false, @islogical);
ip.addParameter('DSR', false, @islogical);
ip.addParameter('xyPixelSize', 0.108, @isnumeric);
ip.addParameter('dz', 0.2, @isnumeric);
ip.addParameter('skewAngle', 32.45, @isnumeric);
ip.addParameter('axisOrder', 'x,y,z', @ischar);
ip.addParameter('dataOrder', 'y,x,z', @ischar);
ip.addParameter('objectiveScan', false, @islogical);
ip.addParameter('IOScan', false, @islogical);
ip.addParameter('zarrFile', false, @islogical);
ip.addParameter('overlapSize', [], @isnumeric);
ip.addParameter('overlapSizeType', 'pixel', @(x) ischar(x) && ismember(lower(x), {'pixel', 'um'}));
ip.addParameter('uuid', '', @ischar);

ip.parse(dataPaths, generationMethod, varargin{:});

pr = ip.Results;
% Overwrite = pr.Overwrite;
channelPatterns = pr.channelPatterns;
tilePatterns = pr.tilePatterns;
xyPixelSize = pr.xyPixelSize;
dz = pr.dz;
skewAngle = pr.skewAngle;
axisOrder = pr.axisOrder;
dataOrder = pr.dataOrder;
DS = pr.DS;
DSR = pr.DSR;
objectiveScan = pr.objectiveScan;
IOScan =  pr.IOScan;
zarrFile = pr.zarrFile;
overlapSize = pr.overlapSize;
overlapSizeType = pr.overlapSizeType;
uuid = pr.uuid;

if ischar(dataPaths)
    dataPaths = {dataPaths};
end
nd = numel(dataPaths);

switch lower(generationMethod)
    case 'encoder'
        for d = 1 : nd
            stitch_generate_imagelist_from_encoder(dataPaths{d}, dz, channelPatterns);
        end
    case 'sqlite'
        for d = 1 : nd        
            stitch_generate_imagelist_from_sqlite(dataPaths{d});        
        end
    case 'tile_position'
        stitch_generate_imagelist_from_tile_positions(dataPaths, channelPatterns=channelPatterns, ...
            tilePatterns=tilePatterns, xyPixelSize=xyPixelSize, dz=dz, skewAngle=skewAngle, ...
            axisOrder=axisOrder, dataOrder=dataOrder, DS=DS, DSR=DSR, objectiveScan=objectiveScan, ...
            IOScan=IOScan, zarrFile=zarrFile, overlapSize=overlapSize, overlapSizeType=overlapSizeType, ...
            uuid=uuid);
    otherwise
        error('Unknown gerneration method, it must be choosen from: encoder, sqlite, or tile_position.');
end

end

