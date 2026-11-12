classdef SettingFilesTest < pktest.TestCase
    % the setting files of an acquisition (parseSettingFile): the stage
    % offset, interval and excitation pixels of XR_parseSettingFiles, and
    % the setting file of each image (XR_parseSettingFiles_wrapper)

    properties (Constant)
        Header = {'***** ***** *****   General   ***** ***** *****', 'Date :  1/1/2024 10:00:00 AM', ''}
        CamName = 'Scan_Iter_%04d_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0000000000msecAbs_000x_000y_000z_0000t'
        NoCamName = 'Scan_Iter_%04d_ch1_CAM1_stack0000_560nm_0000000msec_0000000000msecAbs_000x_000y_000z_0000t'
    end

    methods (Test, TestTags = {'Unit'})
        function lastLineWithoutANewlineIsRead(testCase)
            % the last line read was dropped (feof is true once fgetl has
            % read the last line, with or without a final newline), and
            % the stage line pattern needs a newline after the line: a
            % stage line needed two lines after it
            folder = testCase.createTempFolder();
            fn1 = testCase.writeSettings([folder, 'a_Settings.txt'], [testCase.Header, {testCase.stageLine('-0.5'), 'last line'}], false);
            fn2 = testCase.writeSettings([folder, 'b_Settings.txt'], [testCase.Header, {testCase.stageLine('0.25')}], false);
            fn3 = testCase.writeSettings([folder, 'c_Settings.txt'], [testCase.Header, {testCase.stageLine('0.3'), 'last line'}], true);
            info = XR_parseSettingFiles({fn1, fn2, fn3});
            testCase.verifyEqual([info.StageInterval], [-0.5, 0.25, 0.3]);
            testCase.verifyEqual([info.stageOff], [0, 0, 0]);
            testCase.verifyEqual([info.numExiciation], [101, 101, 101]);
        end

        function onlyTheFirstHundredLinesAreRead(testCase)
            % the stage line can be the 100th line (it was the last line
            % kept, without a newline after it), not the 101st
            folder = testCase.createTempFolder();
            fill = arrayfun(@(x) sprintf('line %d', x), 1 : 96, 'unif', 0);
            fn = testCase.writeSettings([folder, 'a_Settings.txt'], [testCase.Header, fill, {testCase.stageLine('0.5'), 'after'}], true);
            info = XR_parseSettingFiles({fn});
            testCase.verifyEqual(info.StageInterval, 0.5);
            fn = testCase.writeSettings([folder, 'b_Settings.txt'], [testCase.Header, fill, {'one more', testCase.stageLine('0.5'), 'after'}], true);
            testCase.verifyError(@() XR_parseSettingFiles({fn}), ?MException);
        end

        function missingFilesAreNamedInTheWarning(testCase)
            % a missing first file failed (its warning used an undefined
            % name), a missing later file was named as the file before it
            folder = testCase.createTempFolder();
            fn1 = testCase.writeSettings([folder, 'a_Settings.txt'], [testCase.Header, {testCase.stageLine('-0.5'), ''}], true);
            fn3 = testCase.writeSettings([folder, 'c_Settings.txt'], [testCase.Header, {testCase.stageLine('0.5'), ''}], true);
            missing0 = [folder, 'none0_Settings.txt'];
            missing2 = [folder, 'none2_Settings.txt'];

            lastwarn('');
            info = XR_parseSettingFiles({missing0, fn1});
            testCase.verifyEqual(lastwarn(), sprintf('Setting file %s does not exist, skip it.', missing0));
            testCase.verifyEqual(info(2).StageInterval, -0.5);

            lastwarn('');
            info = XR_parseSettingFiles({fn1, missing2, fn3});
            testCase.verifyEqual(lastwarn(), sprintf('Setting file %s does not exist, skip it.', missing2));
            testCase.verifyEqual(info(1).StageInterval, -0.5);
            testCase.verifyEmpty(info(2).StageInterval);
            testCase.verifyEqual(info(3).StageInterval, 0.5);
        end

        function imagesWithAndWithoutACameraInTheirNames(testCase)
            % names with and without a camera were all read with the
            % grammar of names without one: the names with one got the
            % folder as their setting file
            folder = testCase.createTempFolder();
            images = {[folder, sprintf(testCase.CamName, 0), '.tif'], [folder, sprintf(testCase.NoCamName, 1), '.tif'], ...
                [folder, sprintf(testCase.CamName, 2), '.tif'], [folder, sprintf(testCase.NoCamName, 3), '.tif']};
            % by the tile name, by the iteration only, and by the image name
            testCase.writeStageSettings([folder, 'Scan_Iter_0000_000x_000y_000z_0000t_Settings.txt'], '-0.5');
            testCase.writeStageSettings([folder, 'Scan_Iter_0001_000x_000y_000z_0000t_Settings.txt'], '0.5');
            testCase.writeStageSettings([folder, 'Scan_Iter_0002_Settings.txt'], '-0.25');
            testCase.writeStageSettings([folder, sprintf(testCase.NoCamName, 3), '_Settings.txt'], '0.25');
            info = XR_parseSettingFiles_wrapper(images);
            testCase.verifyEqual([info.StageInterval], [-0.5, 0.5, -0.25, 0.25]);
        end

        function imageNameThatMatchesNoGrammar(testCase)
            % a name that matches no grammar got the folder as its setting
            % file (and failed reading it); it has no setting file: a
            % warning, and NaN setting values
            folder = testCase.createTempFolder();
            images = {[folder, sprintf(testCase.CamName, 0), '.tif'], [folder, 'cell_ch0_stack.tif']};
            testCase.writeStageSettings([folder, 'Scan_Iter_0000_000x_000y_000z_0000t_Settings.txt'], '-0.5');
            lastwarn('');
            info = XR_parseSettingFiles_wrapper(images);
            testCase.verifyNumElements(info, 2);
            testCase.verifyEqual(info(1).StageInterval, -0.5);
            testCase.verifyTrue(all(isnan([info(2).stageOff, info(2).StageInterval, info(2).numExiciation])));
            testCase.verifySubstring(lastwarn(), 'cell_ch0_stack.tif');
            % the flags the callers make of them: one per image
            testCase.verifyEqual([info.StageInterval] < 0, [true, false]);

            % as the first image, and alone
            info = XR_parseSettingFiles_wrapper(images([2, 1]));
            testCase.verifyEqual(isnan([info.StageInterval]), [true, false]);
            testCase.verifyEqual(info(2).StageInterval, -0.5);
            info = XR_parseSettingFiles_wrapper(images(2));
            testCase.verifyEqual(fieldnames(info), {'stageOff'; 'StageInterval'; 'numExiciation'});
            testCase.verifyTrue(isnan(info.StageInterval));
        end
    end

    methods (Static, Access = private)
        function [str] = stageLine(interval)
            str = sprintf('X Stage Offset, Interval (um), # of Pixels for Excitation (0) :\t0\t%s\t101', interval);
        end

        function [fn] = writeSettings(fn, lines, finalNewline)
            % a setting file of the lines, with or without a final newline
            str = strjoin(lines, newline);
            if finalNewline
                str = [str, newline];
            end
            fid = fopen(fn, 'w');
            fwrite(fid, str);
            fclose(fid);
        end
    end

    methods (Access = private)
        function [fn] = writeStageSettings(testCase, fn, interval)
            fn = SettingFilesTest.writeSettings(fn, [testCase.Header, {SettingFilesTest.stageLine(interval), ''}], true);
        end
    end
end
