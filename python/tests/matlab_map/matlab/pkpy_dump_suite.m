function pkpy_dump_suite(outFn)
% dump every test of tests/ (also KnownIssue ones) with class, method,
% folder, tags and parameterization, as JSON, for tools/check_test_map.py
% update. Run from the root of an oracle tree after setup([]):
%   matlab -batch "setup([]); addpath('<python>/tests/matlab_map/matlab'); pkpy_dump_suite('suite.json')"
testDir = fullfile(pwd, 'tests');
addpath(testDir);
suite = matlab.unittest.TestSuite.fromFolder(testDir, 'IncludingSubfolders', true);
n = numel(suite);
recs = cell(n, 1);
for i = 1 : n
    t = suite(i);
    ps = t.Parameterization;
    params = struct('property', {}, 'name', {}, 'type', {});
    for j = 1 : numel(ps)
        params(end + 1) = struct('property', ps(j).Property, 'name', ps(j).Name, 'type', class(ps(j))); %#ok<AGROW>
    end
    recs{i} = struct('name', t.Name, 'procedure', t.ProcedureName, 'class', t.TestClass, ...
        'folder', strrep(erase(t.BaseFolder, [testDir, filesep]), '\', '/'), ...
        'tags', {cellstr(t.Tags)}, 'params', params);
end
sel = suite.selectIf(~matlab.unittest.selectors.HasTag('KnownIssue'));
out = struct('matlab', version, 'computer', computer, 'n_all', n, 'n_default', numel(sel), 'tests', {recs});
fid = fopen(outFn, 'w');
fwrite(fid, jsonencode(out, 'PrettyPrint', true));
fclose(fid);
fprintf('n_all %d n_default %d\n', n, numel(sel));
end
