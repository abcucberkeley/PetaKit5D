function [owner] = writeWriterRecord(testCase, resultFn, uuid, owner)
% writes the record of a writer of the result resultFn with uuid, as
% result_writers does (<result folder>/zarr_flag/<name>_<uuid>/writer.txt),
% e.g. of a job of another call that makes the result, dated a minute ago (a
% writer that came first). Without owner (or an empty one), the owner is a
% process that runs on this node during the test (killed at its end), i.e. a
% live writer outside Slurm.

if nargin < 4 || isempty(owner)
    [~, out] = system('sleep 600 > /dev/null 2>&1 & echo $!');
    pid = str2double(strtrim(out));
    testCase.addTeardown(@() system(sprintf('kill %d 2> /dev/null', pid)));
    node = getenv('SLURMD_NODENAME');
    if isempty(node)
        [~, node] = system('uname -n');
        node = strtrim(node);
    end
    owner = sprintf('0_:%s:%d', node, pid);
end
resultFn = strip(regexprep(resultFn, '(?<!^)/+', '/'), 'right', '/');
[resultPath, fsn] = fileparts(resultFn);
flagPath = sprintf('%s/zarr_flag/%s_%s/', resultPath, fsn, uuid);
mkdir(flagPath);
fid = fopen([flagPath, 'writer.txt'], 'w');
fprintf(fid, 'result=%s\nuuid=%s\nowner=%s\n', resultFn, uuid, owner);
fclose(fid);
fileObj = java.io.File([flagPath, 'writer.txt']);
fileObj.setLastModified(fileObj.lastModified() - 60000);

end
