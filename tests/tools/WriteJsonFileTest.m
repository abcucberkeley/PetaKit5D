classdef WriteJsonFileTest < pktest.TestCase
    % writeJsonFile, which writes the parameters.json of every wrapper: it
    % passed the JSON text to fprintf as the format, so a % in a value (e.g.,
    % a file name pattern like x%03d.tif) cut the file off there, and the
    % doubled backslashes of Windows paths and of escaped characters became
    % single ones, which made the file invalid JSON

    methods (Test, TestTags = {'Unit'})
        function fileHoldsTheJsonTextAsItIs(testCase)
            dataPath = testCase.createTempFolder();
            pr = struct('pattern', 'x%03d.tif', 'percent', '100%', 'windowsPath', 'C:\data\new\t1', ...
                'tab', sprintf('a\tb'), 'number', 0.5, 'names', {{'a%s', 'b\n'}});
            fn = [dataPath, 'parameters.json'];
            writeJsonFile(pr, fn);
            testCase.verifyEqual(fileread(fn), jsonencode(pr, PrettyPrint=true));
            read = jsondecode(fileread(fn));
            for f = {'pattern', 'percent', 'windowsPath', 'tab', 'number'}
                testCase.verifyEqual(read.(f{1}), pr.(f{1}), f{1});
            end
            testCase.verifyEqual(reshape(read.names, 1, []), pr.names);
        end
    end
end
