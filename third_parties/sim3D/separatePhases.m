function [] = separatePhases(dataFile, varargin)
% Separate phases of a tiff file into seperate tiff files

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataFile'); % full path to tif file. Ex. /folder/data.tif
ip.addParameter('nphases', 5, @isnumeric);

ip.parse(dataFile, varargin{:});

pr = ip.Results;
nphases = pr.nphases;

data = readtiff(dataFile);
for p=1:nphases
    writetiff(data(:,:,p:nphases:end),[dataFile(1:end-4) '_phase' num2str(p) '.tif']);
end

end
