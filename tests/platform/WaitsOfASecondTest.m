classdef WaitsOfASecondTest < pktest.TestCase
    % pktest.waitsOfASecondProblem, the check of the wrappers' waits between
    % their looks for results (pktest.verifyWaitsOfASecond): the wrapper tests
    % checked that every wait was exactly a second, and failed now and then
    % on slow CI runners (waits of 1 and 1.19 s on GitHub's macOS runner),
    % where wait_for_files waits longer after a slow look, as it should.
    % Whole-round waits between looks, which the wrappers used to make, must
    % still be wrong.

    methods (Test, TestTags = {'Unit'})
        function waitsOfASecondAreRight(testCase)
            testCase.verifyEmpty(pktest.waitsOfASecondProblem([1, 1], 2));
            testCase.verifyEmpty(pktest.waitsOfASecondProblem(ones(1, 7), 7));
            testCase.verifyEmpty(pktest.waitsOfASecondProblem([], 0));
        end

        function waitsAfterSlowLooksAreRight(testCase)
            % ten times a look of 0.1185 s, and the rest of a 3 s round
            testCase.verifyEmpty(pktest.waitsOfASecondProblem([1, 1.185], 2));
            testCase.verifyEmpty(pktest.waitsOfASecondProblem([1, 1.185, 0.815], 3));
        end

        function wholeRoundWaitsAreWrong(testCase)
            testCase.verifySubstring(pktest.waitsOfASecondProblem([3, 3], 2), 'first wait');
            testCase.verifySubstring(pktest.waitsOfASecondProblem([1, 3], 2), 'between 0 and 2 s');
            testCase.verifySubstring(pktest.waitsOfASecondProblem(5, 1), 'first wait');
        end

        function otherCountsAndWaitsAreWrong(testCase)
            testCase.verifySubstring(pktest.waitsOfASecondProblem([1, 1, 1], 2), '3 waits instead of 2');
            testCase.verifySubstring(pktest.waitsOfASecondProblem([], 2), '0 waits instead of 2');
            testCase.verifySubstring(pktest.waitsOfASecondProblem([1.5, 1], 2), 'first wait');
            testCase.verifySubstring(pktest.waitsOfASecondProblem([1, 0], 2), 'between 0 and 2 s');
        end

        function verifyGivesTheWaits(testCase)
            % the wrapper tests' diagnostic names the waits
            testCase.verifySubstring(pktest.waitsOfASecondProblem([1, 1.25], 1), '[1 1.25] s');
        end
    end
end
