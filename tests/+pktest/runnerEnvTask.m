function [] = runnerEnvTask(outFn)
% task function for the task runner tests: writes the metrics environment
% it sees (PETAKIT5D_TASK_METRICS|PETAKIT5D_TASK) to outFn

fid = fopen(outFn, 'w');
fprintf(fid, '%s|%s\n', getenv('PETAKIT5D_TASK_METRICS'), getenv('PETAKIT5D_TASK'));
fclose(fid);

end
