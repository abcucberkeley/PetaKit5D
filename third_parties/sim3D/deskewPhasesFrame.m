function [] = deskewPhasesFrame(dataFile,varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataFile');
ip.addParameter('nphases', 5, @isnumeric);

ip.parse(dataFile, varargin{:});

pr = ip.Results;
nphases = pr.nphases;

separatePhases(dataFile, 'nphases', nphases);

combinePhases();

end

