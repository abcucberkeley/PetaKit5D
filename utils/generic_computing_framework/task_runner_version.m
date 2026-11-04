function [v] = task_runner_version()
% the version of the task runner protocol (task_runner('version'), which
% compile_mccMaster writes with a build), apart so that checking a build
% doesn't load task_runner (36 ms)

v = 1;

end
