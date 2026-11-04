function [logDir, inData] = defaultJobLogDir(dataPath)
% the job log folder the job framework's wrappers use by default (jobLogDir
% '../job_logs', from the PetaKit5D folder, where they submit from):
% ../job_logs beside the PetaKit5D folder if it is there, else
% <dataPath>/job_logs (inData true). As a full path with symbolic links
% resolved, as pktest.fullPath gives it for the paths the jobs get

logDir = fullfile(pktest.repoRoot(), '..', 'job_logs');
inData = ~isfolder(logDir);
if inData
    logDir = fullfile(dataPath, 'job_logs');
end
logDir = pktest.fullPath(logDir);

end
