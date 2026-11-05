classdef UsableGPUNumTest < pktest.TestCase
    % usable_gpu_num, which counts the GPUs PetaKit5D can use and turns on
    % forward compatibility for GPUs newer than this MATLAB supports (e.g.,
    % B200 GPUs, compute capability 10.0, which R2024b doesn't count as
    % available without it, so GPU jobs ran on the CPU). gpuDeviceCount and
    % gpuDevice are stubs that act as such GPUs by the real forward
    % compatibility setting, which the tests turn off and restore, so they
    % run only on machines without GPUs: in R2024b, setting it before CUDA
    % starts can keep turning it on later from taking effect in the session
    % (none of a B200 node's GPUs was then counted, so the GPU tests were
    % skipped). On GPU machines the GPU tests check the real thing.

    methods (Test, TestTags = {'Unit'})
        function gpusNewerThanMATLABGetForwardCompatibility(testCase)
            testCase.stubGPUs(2, 'tooNew');
            [out, n, reason, nAll] = evalc('usable_gpu_num()');
            testCase.verifyEqual([n, nAll], [2, 2]);
            testCase.verifyEmpty(reason);
            testCase.verifyTrue(parallel.gpu.enableCUDAForwardCompatibility());
            % it says so, once: they are available from then on
            testCase.verifySubstring(out, 'forward compatibility is turned on');
            [out, n] = evalc('usable_gpu_num()');
            testCase.verifyEqual(n, 2);
            testCase.verifyEmpty(out);
        end

        function forwardCompatibilityCanBeTurnedOff(testCase)
            % with MATLAB's environment variable at 0
            testCase.stubGPUs(2, 'tooNew');
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('MW_CUDA_FORWARD_COMPATIBILITY', '0'));
            [n, reason, nAll] = usable_gpu_num();
            testCase.verifyEqual([n, nAll], [0, 2]);
            testCase.verifySubstring(reason, 'Selecting GPU 1 failed: GPU device is not supported because it has a higher compute capability');
            testCase.verifyFalse(parallel.gpu.enableCUDAForwardCompatibility());
        end

        function otherFailuresLeaveForwardCompatibilityOff(testCase)
            % e.g., a GPU in use or broken: the reason is MATLAB's, and
            % forward compatibility isn't tried (the GPUs are counted once)
            folder = testCase.stubGPUs(1, 'busy');
            [n, reason, nAll] = usable_gpu_num();
            testCase.verifyEqual([n, nAll], [0, 1]);
            testCase.verifySubstring(reason, 'Selecting GPU 1 failed: All CUDA-capable devices are busy');
            testCase.verifyFalse(parallel.gpu.enableCUDAForwardCompatibility());
            testCase.verifyEqual(strsplit(strtrim(fileread([folder, 'available.txt'])), newline), {'0'});
        end

        function forwardCompatibilityThatDoesNotHelpIsTurnedOffAgain(testCase)
            testCase.stubGPUs(1, 'tooNewBroken');
            [n, reason] = usable_gpu_num();
            testCase.verifyEqual(n, 0);
            testCase.verifySubstring(reason, 'higher compute capability');
            testCase.verifyFalse(parallel.gpu.enableCUDAForwardCompatibility());
        end

        function availableGPUsAreCountedAsTheyAre(testCase)
            folder = testCase.stubGPUs(3, 'available');
            [n, reason, nAll] = usable_gpu_num();
            testCase.verifyEqual([n, nAll], [3, 3]);
            testCase.verifyEmpty(reason);
            testCase.verifyFalse(parallel.gpu.enableCUDAForwardCompatibility());
            testCase.verifyFalse(isfile([folder, 'selected.txt']), 'A GPU was selected.');
        end

        function noGPUFound(testCase)
            folder = testCase.stubGPUs(0, 'available');
            [n, reason, nAll] = usable_gpu_num();
            testCase.verifyEqual([n, nAll], [0, 0]);
            testCase.verifyEmpty(reason);
            testCase.verifyFalse(isfile([folder, 'selected.txt']), 'A GPU was selected.');
        end
    end

    methods (Access = private)
        function [folder] = stubGPUs(testCase, found, kind)
            % found GPUs of a kind: 'available'; 'tooNew' (available only
            % with forward compatibility; selecting one without it fails with
            % MATLAB's DeviceTooNew error); 'tooNewBroken' (as 'tooNew', but
            % not available with it either); 'busy' (not available, and
            % selecting one fails otherwise). In the returned folder,
            % gpuDevice writes selected.txt at each call, and gpuDeviceCount
            % adds the forward compatibility setting to available.txt at each
            % count of the available GPUs.
            nReal = 0;
            try
                nReal = gpuDeviceCount("all");
            catch
            end
            testCase.assumeEqual(nReal, 0, 'The stubs stand in for GPUs on machines without them.');
            try
                fc = parallel.gpu.enableCUDAForwardCompatibility();
            catch ME
                testCase.assumeFail(sprintf('The forward compatibility setting is not available here: %s', ME.message));
            end
            testCase.assumeFalse(fc, 'Forward compatibility is on in this MATLAB session.');
            testCase.addTeardown(@() parallel.gpu.enableCUDAForwardCompatibility(fc));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('MW_CUDA_FORWARD_COMPATIBILITY', ''));

            folder = testCase.createTempFolder();
            fid = fopen([folder, 'gpuDeviceCount.m'], 'w');
            fprintf(fid, '%s\n', 'function [n] = gpuDeviceCount(varargin)', ...
                '% test stub written by UsableGPUNumTest', sprintf('n = %d;', found), ...
                'if nargin > 0 && strcmp(varargin{1}, ''available'')', sprintf('    kind = ''%s'';', kind), ...
                '    fc = parallel.gpu.enableCUDAForwardCompatibility();', ...
                '    fid = fopen([fileparts(mfilename(''fullpath'')), ''/available.txt''], ''a'');', ...
                '    fprintf(fid, ''%d\n'', fc);', '    fclose(fid);', ...
                '    n = n * (strcmp(kind, ''available'') || (strcmp(kind, ''tooNew'') && fc));', ...
                'end', 'end');
            fclose(fid);
            fid = fopen([folder, 'gpuDevice.m'], 'w');
            fprintf(fid, '%s\n', 'function [varargout] = gpuDevice(varargin)', ...
                '% test stub written by UsableGPUNumTest: selecting a GPU fails as for the kind of GPU', ...
                sprintf('kind = ''%s'';', kind), ...
                'fid = fopen([fileparts(mfilename(''fullpath'')), ''/selected.txt''], ''a'');', 'fprintf(fid, ''1\n'');', 'fclose(fid);', ...
                'if any(strcmp(kind, {''tooNew'', ''tooNewBroken''})) && ~parallel.gpu.enableCUDAForwardCompatibility()', ...
                ['    error(''parallel:gpu:device:DeviceTooNew'', ''GPU device is not supported because it has a higher ' ...
                'compute capability (10.0) than is natively supported by the CUDA libraries (test stub).'');'], ...
                'elseif ~strcmp(kind, ''available'')', ...
                '    error(''parallel:gpu:device:DeviceBusy'', ''All CUDA-capable devices are busy or unavailable (test stub).'');', ...
                'end', 'varargout = {struct(''Index'', 1)};', 'end');
            fclose(fid);
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture('MATLAB:dispatcher:nameConflict'));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(folder));
            clear gpuDeviceCount gpuDevice
            testCase.addTeardown(@() clear('gpuDeviceCount', 'gpuDevice'));
        end
    end
end
