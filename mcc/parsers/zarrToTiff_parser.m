function [] = zarrToTiff_parser(zarrFilename, tiffFilename, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('tifFilename', @ischar);
ip.addRequired('zarrFilename', @ischar);
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('usrFcn', '', @(x) isempty(x) || isa(x,'function_handle') || ischar(x));
ip.addParameter('uuid', '', @ischar);

ip.parse(zarrFilename, tiffFilename, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
usrFcn = pr.usrFcn;
uuid = pr.uuid;

if ischar(Overwrite)
    Overwrite = str2num(Overwrite);
end

zarrToTiff(tifFilename, zarrFilename, Overwrite=Overwrite, usrFcn=usrFcn, uuid=uuid);

end

