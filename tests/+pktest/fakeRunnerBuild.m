function [b] = fakeRunnerBuild(folder, withRunners)
% a stand-in compiled build for task runner tests: run_mccMaster.sh as the
% compiler makes it (its argument handling, Linux) and a bash mccMaster that
% runs the task 'pk_task <output> <mode>' (mode ok: writes x to the output;
% slow: the same after 2 s; error: fails; crash: kills its process)
% and, as a task runner (PETAKIT5D_RUNNER set), claims and runs the next
% tasks through the job script as task_runner does. The build has task
% runners of this code's version (task_runner.txt) unless withRunners is
% false. Fields: master (run_mccMaster.sh), mcr (a MATLAB Runtime folder),
% log (each task run: output, mode, PETAKIT5D_TASK or -, process, CUDA
% device, runner 0/1).

if nargin < 2
    withRunners = true;
end
if ~endsWith(folder, '/')
    folder = [folder, '/'];
end
mkdir(folder);
b.master = [folder, 'run_mccMaster.sh'];
b.mcr = [folder, 'mcr'];
b.log = [folder, 'runs.log'];
mkdir([b.mcr, '/bin']);
write_lines(b.master, {
    '#!/bin/sh'
    'exe_name=$0'
    'exe_dir=`dirname "$0"`'
    'if [ "x$1" = "x" ]; then'
    '  echo Usage:'
    'else'
    '  MCRROOT="$1"'
    '  LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;'
    '  export LD_LIBRARY_PATH;'
    '  shift 1'
    '  args='
    '  while [ $# -gt 0 ]; do'
    '      token=$1'
    '      args="${args} \"${token}\"" '
    '      shift'
    '  done'
    '  eval "\"${exe_dir}/mccMaster\"" $args'
    'fi'
    'exit'});
write_lines([folder, 'mccMaster'], {
    '#!/bin/bash'
    'CFG="$PETAKIT5D_RUNNER"'
    'R=0; [ -n "$CFG" ] && R=1'
    'unset PETAKIT5D_RUNNER'
    'run_one() {'
    sprintf('    echo "$2 $3 ${PETAKIT5D_TASK:--} $$ ${CUDA_VISIBLE_DEVICES:-none} $R" >> "%s"', b.log)
    '    case "$3" in'
    '        error) return 1 ;;'
    '        crash) kill -9 $$ ;;'
    '        slow) sleep 2; echo x > "$2" ;;'
    '        *) echo x > "$2" ;;'
    '    esac'
    '}'
    'run_one "$@" || exit 1'
    '[ -z "$CFG" ] && exit 0'
    'IFS="|" read -r D own steal s st L task <<< "$CFG"'
    'printf "d %s\nt 0.1\n" "$task" >> "$st"'
    'while true; do'
    '    bash "$D/run_tasks.sh" claim "$own" "$steal" "$s" 1 "$st" "$L" "$L/next_$s"'
    '    f=$(cat "$L/next_$s")'
    '    [ -z "$f" ] && { echo "x none" >> "$st"; exit 0; }'
    '    echo "c $f" >> "$st"'
    '    mapfile -d "" parts < "$L/dump_$s"'
    '    if PETAKIT5D_TASK="${parts[1]}" run_one "${parts[@]:3}"; then'
    '        printf "d %s\nt 0.1\n" "$f" >> "$st"'
    '    else'
    '        printf "e %s\nx failed\n" "$f" >> "$st"'
    '        exit 0'
    '    fi'
    'done'});
fileattrib(b.master, '+x');
fileattrib([folder, 'mccMaster'], '+x');
if withRunners
    write_lines([folder, 'task_runner.txt'], {sprintf('%d', task_runner('version'))});
end

end


function [] = write_lines(fn, lines)

fid = fopen(fn, 'w');
fwrite(fid, [strjoin(lines, newline), newline]);
fclose(fid);

end
