function [] = XR_resample_dataset_parser(dataPath, resultPath, rsfactor, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPath', @ischar);
ip.addRequired('resultPath', @ischar);
ip.addRequired('rsfactor', @(x) isnumeric(x) || ischar(x));
ip.addParameter('Interp', 'linear', @ischar);
ip.addParameter('Save16bit', true, @(x) islogical(x) || ischar(x));
ip.addParameter('zarrFile', false, @(x) islogical(x) || ischar(x));
ip.addParameter('saveZarr', false, @(x) islogical(x) || ischar(x)); % use zarr file as output
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x)); % master node participate in the task computing. 
ip.addParameter('cpusPerTask', 2, @(x) isnumeric(x) || ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ConfigFile', '', @ischar);

ip.parse(dataPath, resultPath, rsfactor, varargin{:});

pr = ip.Results;
Interp = pr.Interp;
Save16bit = pr.Save16bit;
zarrFile = pr.zarrFile;
saveZarr = pr.saveZarr;
parseCluster = pr.parseCluster;
jobLogDir = pr.jobLogDir;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

if ischar(rsfactor)
    rsfactor = str2num(rsfactor);
end
if ischar(Save16bit)
    Save16bit = str2num(Save16bit);
end
if ischar(zarrFile)
    zarrFile = str2num(zarrFile);
end
if ischar(saveZarr)
    saveZarr = str2num(saveZarr);
end
if ischar(parseCluster)
    parseCluster = str2num(parseCluster);
end
if ischar(masterCompute)
    masterCompute = str2num(masterCompute);
end
if ischar(cpusPerTask)
    cpusPerTask = str2num(cpusPerTask);
end
if ischar(mccMode)
    mccMode = str2num(mccMode);
end

XR_resample_dataset(dataPath, resultPath, rsfactor, Interp=Interp, Save16bit=Save16bit, ...
    zarrFile=zarrFile, saveZarr=saveZarr, parseCluster=parseCluster, jobLogDir=jobLogDir, ...
    masterCompute=masterCompute, cpusPerTask=cpusPerTask, uuid=uuid, mccMode=mccMode, ...
    ConfigFile=ConfigFile);

end

