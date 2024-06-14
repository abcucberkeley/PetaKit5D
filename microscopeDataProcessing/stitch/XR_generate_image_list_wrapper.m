function [] = XR_generate_image_list_wrapper(dataPaths, generationMethod, varargin)
% wrapper to generate image list for a given generation method: encoder,
% sqlite, or tile_position
% 
% It saves as a csv file with name ImageList_from_*.csv in the dataPath. The format is consistent with old csv files. 
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
ip.addParameter('xyPixelSize', 0.108, @isnumeric); % in um
ip.addParameter('dz', 0.2, @isnumeric); % in um
ip.addParameter('skewAngle', 32.45, @isnumeric);
ip.addParameter('axisOrder', 'x,y,z', @ischar); % tile indices axis order, -x means scan in a reversed direction in x axis
ip.addParameter('dataOrder', 'y,x,z', @ischar); % data axis order
ip.addParameter('objectiveScan', false, @islogical);
ip.addParameter('IOScan', false, @islogical);
ip.addParameter('zarrFile', false, @islogical);
ip.addParameter('overlapSize', [], @isnumeric); % in stage coordinates
ip.addParameter('overlapSizeType', 'pixel', @(x) ischar(x) && ismember(lower(x), {'pixel', 'um'})); % yxz in pixel, xyz in um
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

