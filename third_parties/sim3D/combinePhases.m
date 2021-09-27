function [] = combinePhases(dataFile, varargin)
% Combine phases into one tiff file

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataFile'); % full path to combined phases without phase end tag. Ex. /folder/data.tif combines files like /folder/data_phase1.tif
ip.addParameter('nphases', 5, @isnumeric);

ip.parse(dataFile, varargin{:});

pr = ip.Results;
nphases = pr.nphases;

% get first phase to check size
cPhase = readtiff([dataFile(1:end-4) '_phase' num2str(1) '.tif']);
outSize = size(cPhase);

out = zeros(outSize(1),outSize(2),outSize(3)*nphases);
for p=1:nphases
    if(p > 1)
        cPhase = readtiff([dataFile(1:end-4) '_phase' num2str(p) '.tif']);
    end
    out(:,:,p:nphases:end) = cPhase;
end

writetiff(single(out),[dataFile(1:end-4) '_combined.tif']);

end