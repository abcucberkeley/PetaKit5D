classdef MccBuildCheckTest < pktest.TestCase
    % compile_mccMaster saves the task interfaces of a build (the argument
    % names of mcc/parsers/*_parser.m) in parsers_manifest.txt, and
    % check_mcc_master_build compares them with the code before compiled
    % tasks are submitted: a build made before a task got a new option
    % rejects it, so every task of a run failed (the reason only in the job
    % logs), e.g., after 08e857d or f667138 with an older build.

    methods (Test, TestTags = {'Unit'})
        function manifestListsTheArgumentsOfEachTask(testCase)
            manifest = splitlines(strtrim(mcc_parser_manifest(testCase.parserDir())));
            parsers = dir([testCase.parserDir(), '/*_parser.m']);
            testCase.verifyNumElements(manifest, numel(parsers));
            line = manifest(startsWith(manifest, 'XR_crop_block: '));
            testCase.verifyEqual(line, {['XR_crop_block: batchInds,zarrFullpath,cropFullpath,flagFullname,', ...
                'BatchBBoxes,RegionBBoxes,Overwrite,mipBatchSize,uuid,debug']});
        end

        function buildOfThisCodePasses(testCase)
            masterFn = testCase.writeBuild(mcc_parser_manifest(testCase.parserDir()));
            testCase.verifyWarningFree(@() check_mcc_master_build(masterFn));
        end

        function buildWithOtherTaskInterfacesStops(testCase)
            % a build from before XR_crop_block got mipBatchSize
            manifest = strrep(mcc_parser_manifest(testCase.parserDir()), 'Overwrite,mipBatchSize,uuid', 'Overwrite,uuid');
            masterFn = testCase.writeBuild(manifest);
            testCase.verifyError(@() check_mcc_master_build(masterFn), 'mcc:outdatedBuild');
            try
                check_mcc_master_build(masterFn);
            catch ME
                testCase.verifySubstring(ME.message, 'XR_crop_block');
                testCase.verifySubstring(ME.message, 'compile_mccMaster');
            end
        end

        function buildWithoutManifestWarns(testCase)
            % builds made before the check still run
            masterFn = testCase.writeBuild('');
            delete([fileparts(masterFn), '/parsers_manifest.txt']);
            testCase.verifyWarning(@() check_mcc_master_build(masterFn), 'mcc:unverifiedBuild');
        end
    end

    methods (Access = private)
        function [masterFn] = writeBuild(testCase, manifest)
            % a build folder with run_mccMaster.sh and the given manifest
            buildDir = testCase.createTempFolder();
            masterFn = [buildDir, 'run_mccMaster.sh'];
            fclose(fopen(masterFn, 'w'));
            fid = fopen([buildDir, 'parsers_manifest.txt'], 'w');
            fprintf(fid, '%s', manifest);
            fclose(fid);
        end
    end

    methods (Static, Access = private)
        function [d] = parserDir()
            d = fileparts(which('mccMaster'));
        end
    end
end
