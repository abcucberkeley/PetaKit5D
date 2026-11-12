function [] = verifyWaitsOfASecond(testCase, pauses, n, msg)
% verifies the waits of a wrapper between its looks for results or inputs,
% as the pause calls a FunctionSpyFixture('pause') recorded: n waits, the
% first of a second and none of 2 s or more (pktest.waitsOfASecondProblem
% says why the waits are not all exactly a second).
%
% msg: the diagnostic of the check (default: none).


if nargin < 4
    msg = '';
end
waits = cellfun(@(x) x{1}, pauses.getCalls());
problem = pktest.waitsOfASecondProblem(waits, n);
testCase.verifyEmpty(problem, strtrim(sprintf('%s %s', msg, problem)));

end
