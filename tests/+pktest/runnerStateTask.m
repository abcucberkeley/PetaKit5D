function [] = runnerStateTask(outFn, mode)
% task function for the task runner tests: saves the process state it sees
% to outFn (a .mat file), then changes it all the ways a task can, so the
% next task in the same process shows whether the runner reset it:
% current folder, path, warning states and last warning, random numbers,
% environment, global and base workspace variables, a persistent variable,
% open files, figures. mode 'error' throws after saving; 'hold' also keeps
% memory in root appdata (not reset, for the memory guard); 'threads' also
% changes the number of computational threads (the process then ends).

if nargin < 2
    mode = 'ok';
end

persistent nCalls
if isempty(nCalls)
    nCalls = 0;
end
nCalls = nCalls + 1;

s.pwd = pwd;
s.pathLength = numel(path);
s.warningState = warning('query', 'pktest:runnerState');
s.warningState = s.warningState.state;
[s.lastMessage, s.lastId] = lastwarn;
s.rand = rand(1, 3);
s.threads = maxNumCompThreads;
s.env = getenv('PK_RUNNER_STATE');
s.runnerEnv = getenv('PETAKIT5D_RUNNER');
s.hasGlobal = ~isempty(who('global', 'pkRunnerStateGlobal'));
s.hasBase = evalin('base', 'exist(''pkRunnerStateBase'', ''var'')');
s.persistentCalls = nCalls;
s.openFiles = numel(fopen_all());
s.figures = 0;
try
    s.figures = numel(groot().Children);
catch
end
save(outFn, '-struct', 's');

if strcmp(mode, 'error')
    error('pktest:runnerStateTask:failed', 'Task failed on purpose for %s', outFn);
end

% changes
d = tempname;
mkdir(d);
cd(d);
addpath(d);
warning('off', 'pktest:runnerState');
w = warning('off', 'backtrace');
warning('pktest:runnerStateLast', 'a warning of the task');
warning(w);
rand(1, 100);
if strcmp(mode, 'threads')
    maxNumCompThreads(1 + (maxNumCompThreads == 1));
end
setenv('PK_RUNNER_STATE', 'changed');
global pkRunnerStateGlobal %#ok<GVMIS>
pkRunnerStateGlobal = 1;
assignin('base', 'pkRunnerStateBase', 1);
fopen([d, '/open.txt'], 'w');
try
    % (no figures without Java)
    figure('Visible', 'off');
catch
end
if strcmp(mode, 'hold')
    % 2 GB (well above the test's margin, as memory other tests left can be
    % freed at the same time)
    setappdata(0, 'pkRunnerStateHold', ones(1, 2^28));
end

end


function [fids] = fopen_all()

try
    fids = openedFiles();
catch
    fids = fopen('all');
end

end
