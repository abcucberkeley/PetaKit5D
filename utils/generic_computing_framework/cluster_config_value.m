function [value] = cluster_config_value(configFile, name, value)
% the value of the parameter name in the cluster configuration file
% configFile (json or mat, as the job framework reads it; keys match names
% regardless of case), or the given value if there is no file or it has no
% such key

if isempty(configFile) || ~exist(configFile, 'file')
    return;
end
[~, ~, ext] = fileparts(configFile);
switch ext
    case '.json'
        confData = jsondecode(fileread(configFile));
    case '.mat'
        confData = load(configFile);
    otherwise
        return;
end
keys = fieldnames(confData);
key = keys(strcmpi(keys, name));
if ~isempty(key)
    value = confData.(key{end});
end

end
