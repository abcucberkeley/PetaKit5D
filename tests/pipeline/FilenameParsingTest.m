classdef FilenameParsingTest < pktest.TestCase
    % the file lists of the pipeline: channel patterns, partial volume files
    % (groupPartialVolumeFiles, XR_parseImageFilenames) and the grouping of
    % the files of a time point by channel (groupImageFilenamesByChannels)

    methods (Test, TestTags = {'Unit'})
        function partsOfNamesWithRegexpCharacters(testCase)
            % the main names were used as regular expressions, so the parts
            % of a name with '+' or '(1)' were lost
            dataPath = testCase.createTempFolder();
            names = {'cell(1)_CamA', 'cell(1)_CamA_part0001', 'y+z_CamA', 'y+z_CamA_part0001', 'y+z_CamA_part0002'};
            depths = [3, 2, 4, 2, 1];
            for i = 1 : numel(names)
                writetiff(zeros(8, 6, depths(i), 'uint16'), [dataPath, names{i}, '.tif']);
            end
            [fnames, ~, gfnames, partialvols, dataSizes] = XR_parseImageFilenames({dataPath}, {'CamA'}, ...
                false, false, false, {}, false, 1, false);
            testCase.verifyEqual(fnames, {'cell(1)_CamA.tif'; 'y+z_CamA.tif'});
            testCase.verifyEqual(gfnames, {{'cell(1)_CamA.tif'; 'cell(1)_CamA_part0001.tif'}; ...
                {'y+z_CamA.tif'; 'y+z_CamA_part0001.tif'; 'y+z_CamA_part0002.tif'}});
            testCase.verifyEqual(partialvols, [true; true]);
            testCase.verifyNumElements(dataSizes, 2);
        end

        function partialVolumeFlagsAreLogicalWithAnEmptyFolder(testCase)
            % partialvols was double when a folder had no files
            emptyPath = testCase.createTempFolder();
            dataPath = testCase.createTempFolder();
            writetiff(zeros(8, 6, 3, 'uint16'), [dataPath, 'im_CamA.tif']);
            [fnames, fdinds, ~, partialvols] = XR_parseImageFilenames({emptyPath, dataPath}, {'CamA'}, ...
                false, false, false, {}, false, 1, false);
            testCase.verifyEqual(fnames, {'im_CamA.tif'});
            testCase.verifyEqual(fdinds, 2);
            testCase.verifyEqual(partialvols, false);
        end

        function fileListAsARow(testCase)
            % a file list given as a row failed (rows of part names put
            % under the main name)
            dataPath = testCase.createTempFolder();
            names = {'a_CamA.tif', 'a_CamA_part0002.tif', 'a_CamA_part0001.tif', 'b_CamA.tif'};
            for i = 1 : numel(names)
                fclose(fopen([dataPath, names{i}], 'w'));
            end
            fns = strcat(dataPath, names);
            [containPart, groups] = groupPartialVolumeFiles('', fns);
            [containPartCol, groupsCol, datenums, datasizes] = groupPartialVolumeFiles('', fns(:));
            expected = {{'a_CamA.tif'; 'a_CamA_part0001.tif'; 'a_CamA_part0002.tif'}; {'b_CamA.tif'}};
            testCase.verifyEqual(groups, expected);
            testCase.verifyEqual(groupsCol, expected);
            testCase.verifyEqual(containPart(:), containPartCol);
            testCase.verifyEqual(cellfun(@numel, datenums), [3; 1]);
            testCase.verifyEqual(cellfun(@numel, datasizes), [3; 1]);
        end

        function channelPatternThatIsNotARegexp(testCase)
            % contains(names, regexpPattern(p)) fails for a pattern with an
            % unmatched '(', '[' or '{', also where the names contain it
            dataPath = testCase.createTempFolder();
            names = {'im_Cam(A_ch0.tif', 'im_Cam[B_ch0.tif', 'im_CamC_ch0.tif'};
            for i = 1 : numel(names)
                writetiff(zeros(8, 6, 3, 'uint16'), [dataPath, names{i}]);
            end
            fnames = XR_parseImageFilenames({dataPath}, {'Cam(A', 'Cam[B'}, false, false, false, {}, false, 1, false);
            testCase.verifyEqual(fnames, names(1 : 2)');

            [fnames, ~, ~, ~, chInds] = parseImageFilenames(dataPath, false, {'Cam[B', 'cam(a'});
            testCase.verifyEqual(fnames, names(1 : 2)');
            testCase.verifyEqual(chInds, [2; 1]);

            testCase.verifyEqual(groupImageFilenamesByChannels(names(1 : 2), {'Cam(A', 'Cam[B'}), [1, 2]);
        end

        function channelPatternsStayRegularExpressions(testCase)
            % valid patterns are still regular expressions, others text
            testCase.verifyEqual(contains({'a_ch1_b', 'a_chX_b'}, channel_regexp_pattern('_ch[0-9]_')), [true, false]);
            testCase.verifyEqual(contains({'a_Cam{A', 'a_CamA'}, channel_regexp_pattern('Cam{A')), [true, false]);
        end

        function channelsWithDifferentNumbersOfFiles(testCase)
            % a channel without a file left for a group failed (assignment of
            % an empty max); the group gets 0 for it
            fsn = 'Scan_Iter_0000_Cam%s_ch0_t%04d.tif';
            names = {sprintf(fsn, 'A', 1), sprintf(fsn, 'A', 2), sprintf(fsn, 'A', 3), ...
                sprintf(fsn, 'B', 1), sprintf(fsn, 'B', 2)};
            testCase.verifyEqual(groupImageFilenamesByChannels(names, {'CamA', 'CamB'}, [1, 1, 1, 2, 2]), ...
                [1, 4; 2, 5; 3, 0]);
            testCase.verifyEqual(groupImageFilenamesByChannels(names([4, 1, 5, 2]), {'CamA', 'CamB'}, [2, 1, 2, 1]), ...
                [2, 1; 4, 3]);
            % more files of the second channel
            names = names([1, 2, 4, 5]);
            names{end + 1} = sprintf(fsn, 'B', 3);
            testCase.verifyEqual(groupImageFilenamesByChannels(names, {'CamA', 'CamB'}), [1, 3; 2, 4; 0, 5]);
            % the same numbers: as before
            testCase.verifyEqual(groupImageFilenamesByChannels(names(1 : 4), {'CamA', 'CamB'}), [1, 3; 2, 4]);
        end
    end
end
