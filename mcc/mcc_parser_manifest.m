function [manifest] = mcc_parser_manifest(parserDir)
% the task interfaces of the compiled master: one line per parser file in
% parserDir (<function>_parser.m) with the names of its required, optional
% and name-value arguments in order, e.g., 'XR_crop_block: batchInds,...,uuid,debug'.
% compile_mccMaster saves it with a build (parsers_manifest.txt), and
% check_mcc_master_build compares the build's with the code's.

dir_info = dir([parserDir, '/*_parser.m']);
fns = sort({dir_info.name});
lines = cell(numel(fns), 1);
for i = 1 : numel(fns)
    txt = fileread([parserDir, '/', fns{i}]);
    % without commented-out lines
    txt = regexprep(txt, '(?m)^\s*%.*$', '');
    args = regexp(txt, 'ip\.add(?:Required|Optional|Parameter)\(\s*''(\w+)''', 'tokens');
    args = cellfun(@(x) x{1}, args, 'unif', 0);
    lines{i} = sprintf('%s: %s', fns{i}(1 : end - numel('_parser.m')), strjoin(args, ','));
end
manifest = sprintf('%s\n', lines{:});

end
