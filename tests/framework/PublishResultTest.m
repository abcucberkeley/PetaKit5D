classdef PublishResultTest < pktest.TestCase
    % keep_other_result and publish_result: a task that is about to move its
    % temporary result into place keeps a result another writer of the same
    % result already put there (two job sets on the same data), instead of
    % deleting it while readers may use it

    methods (Test, TestTags = {'Unit'})
        function resultIsMovedIntoPlace(testCase)
            [tmpFn, resultFn, flagPath] = testCase.makeResult(false);
            testCase.verifyTrue(publish_result(tmpFn, resultFn, flagPath));
            testCase.verifyEqual(testCase.listing(resultFn), {'mine.txt'});
            testCase.verifyFalse(isfolder(tmpFn));
            testCase.verifyFalse(isfolder(flagPath));
        end

        function resultOfAnotherWriterIsKept(testCase)
            [tmpFn, resultFn, flagPath] = testCase.makeResult(true);
            published = [];
            out = evalc('published = publish_result(tmpFn, resultFn, flagPath);');
            testCase.verifyFalse(published);
            testCase.verifyEqual(testCase.listing(resultFn), {'other.txt'}, 'The other writer''s result was replaced.');
            testCase.verifyFalse(isfolder(tmpFn), 'The temporary result was not removed.');
            testCase.verifyFalse(isfolder(flagPath), 'The flags were not removed.');
            testCase.verifySubstring(out, 'made by another job');
        end

        function resultOfAWriterJustBeforeIsKept(testCase)
            % a result that comes between the check and the move: movefile
            % puts the temporary result into it, and it is removed there.
            % The stub makes the check miss the other result.
            [tmpFn, resultFn, flagPath] = testCase.makeResult(true);
            testCase.applyFixture(pktest.FunctionStubFixture('batch_file_exist', false));
            published = [];
            evalc('published = publish_result(tmpFn, resultFn, flagPath);');
            testCase.verifyFalse(published);
            testCase.verifyEqual(testCase.listing(resultFn), {'other.txt'});
            testCase.verifyFalse(isfolder(tmpFn));
            testCase.verifyFalse(isfolder(flagPath));
        end

        function keptResultRemovesATemporaryFile(testCase)
            folder = testCase.createTempFolder();
            resultFn = [folder, 'result.tif'];
            tmpFn = [folder, 'result_1234.tif'];
            fclose(fopen(resultFn, 'w'));
            fclose(fopen(tmpFn, 'w'));
            kept = [];
            evalc('kept = keep_other_result(resultFn, tmpFn);');
            testCase.verifyTrue(kept);
            testCase.verifyFalse(isfile(tmpFn));
            testCase.verifyTrue(isfile(resultFn));
            testCase.verifyFalse(keep_other_result([folder, 'none.tif'], [folder, 'none_1234.tif']));
        end
    end

    methods (Access = private)
        function [tmpFn, resultFn, flagPath] = makeResult(testCase, otherWriter)
            % a temporary result (a folder with mine.txt) and its flags, and
            % with otherWriter the result of another writer (other.txt)
            folder = testCase.createTempFolder();
            resultFn = [folder, 'result.zarr'];
            tmpFn = [resultFn, '_1234'];
            flagPath = [folder, 'zarr_flag/result_1234/'];
            mkdir(tmpFn);
            fclose(fopen([tmpFn, '/mine.txt'], 'w'));
            mkdir(flagPath);
            fclose(fopen([flagPath, 'blocks_1_1.mat'], 'w'));
            if otherWriter
                mkdir(resultFn);
                fclose(fopen([resultFn, '/other.txt'], 'w'));
            end
        end

        function [names] = listing(~, folder)
            d = dir(folder);
            names = setdiff({d.name}, {'.', '..'});
        end
    end
end
