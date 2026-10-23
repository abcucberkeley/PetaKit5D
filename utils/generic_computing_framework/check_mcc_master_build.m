function [] = check_mcc_master_build(MCCMasterStr)
% stop with an error if the compiled master MCCMasterStr was built from code
% with other task interfaces than this code (the mcc parsers): its tasks
% would reject options this code gives them, and each of them would fail.
% compile_mccMaster saves the interfaces of a build in parsers_manifest.txt
% next to run_mccMaster.sh; a build without it (made before the check) only
% gets a warning. Nothing is checked when this code runs compiled itself
% (its tasks come from the same build). A build is checked once per session.

persistent checkedBuilds
if isdeployed
    return;
end
if isempty(checkedBuilds)
    checkedBuilds = {};
end

manifestFn = [fileparts(MCCMasterStr), '/parsers_manifest.txt'];
dir_info = dir(manifestFn);
key = MCCMasterStr;
if ~isempty(dir_info)
    key = sprintf('%s|%.10f', MCCMasterStr, dir_info.datenum);
end
if any(strcmp(checkedBuilds, key))
    return;
end

parserDir = fileparts(which('mccMaster'));
if isempty(dir_info) || isempty(parserDir)
    warning('mcc:unverifiedBuild', ['The compiled master %s has no parsers_manifest.txt (it was built before compile_mccMaster ', ...
        'saved one), so it can''t be checked against this code. If its tasks fail, rebuild it with compile_mccMaster.'], ...
        MCCMasterStr);
else
    buildLines = splitlines(strtrim(fileread(manifestFn)));
    codeLines = splitlines(strtrim(mcc_parser_manifest(parserDir)));
    if ~isequal(buildLines, codeLines)
        funcs = unique(regexprep(setxor(buildLines, codeLines), ':.*$', ''), 'stable');
        error('mcc:outdatedBuild', ['The compiled master %s is out of date: the task interfaces of %s ', ...
            'differ from this code. Rebuild it with compile_mccMaster.'], MCCMasterStr, ...
            strjoin(funcs(1 : min(3, end)), ', '));
    end
end
checkedBuilds{end + 1} = key;

end
