classdef PSFAnalysisTest < pktest.TestCase
    % XR_psf_analysis_wrapper ("PSF analysis" in demo_useful_tools.m):
    % deskew the PSF, then plot it and its OTF against a Richards & Wolf
    % widefield PSF. Even for a 64 x 64 x 101 PSF this needs ~22 GB of RAM
    % (the plots are simulated on a 1001^3 double grid) and ~80 s, so it is
    % skipped on machines with less than 32 GB (e.g., GitHub runners).

    methods (Test, TestTags = {'Integration', 'Slow'})
        function analysisWritesFiguresAndInfo(testCase)
            testCase.assumeGreaterThanOrEqual(testCase.memoryGB(), 32, 'PSF analysis needs ~22 GB of RAM.');
            dataPath = testCase.createTempFolder();
            psfPath = [dataPath, 'PSF/'];
            mkdir([psfPath, 'RW_PSFs/']);
            fsn = '488_CamA_ch0';
            writetiff(uint16(pktest.beads([64, 64, 101], [32, 32, 51], [1.6, 3], 2000, 100)), [psfPath, fsn, '.tif']);
            rwFn = [psfPath, 'RW_PSFs/rw.tif'];
            writetiff(single(pktest.beads([128, 128, 101], [64, 64, 51], [2, 6], 1, 0)), rwFn);

            XR_psf_analysis_wrapper({psfPath}, 'xyPixelSize', 0.108, 'dz', 0.3, 'skewAngle', 32.45, ...
                'channelPatterns', {'CamA_ch0'}, 'channels', 488, 'deskew', true, 'RWFn', {rwFn}, ...
                'sourceStr', 'test', 'visible', false, 'parseCluster', false, 'masterCompute', true);

            resultPath = [psfPath, 'DS/PSFAnalysis/'];
            % deskewed PSF (the demo's PSF is in skewed space)
            ds = readtiff([psfPath, 'DS/', fsn, '.tif']);
            testCase.verifyEqual(size(ds, [1, 3]), [64, 101]);
            testCase.verifyGreaterThan(size(ds, 2), 64);
            for f = {['wT_', fsn, '.png'], ['comp_', fsn, '.png'], [fsn, '_infos.mat'], 'RW_info.mat'}
                testCase.verifyTrue(exist([resultPath, f{1}], 'file') == 2, sprintf('Missing %s.', f{1}));
            end
        end
    end

    methods (Access = private)
        function [gb] = memoryGB(~)
            % memory available to this MATLAB: the Slurm allocation if any,
            % otherwise the machine's memory
            slurmMB = str2double(getenv('SLURM_MEM_PER_NODE'));
            if ~isnan(slurmMB)
                gb = slurmMB / 1024;
            elseif ispc
                m = memory;
                gb = m.MemAvailableAllArrays / 1024 ^ 3;
            elseif exist('/proc/meminfo', 'file')
                kb = regexp(fileread('/proc/meminfo'), 'MemTotal:\s+(\d+)', 'tokens', 'once');
                gb = str2double(kb{1}) / 1024 ^ 2;
            else
                gb = getSystemMemory();
            end
        end
    end
end
