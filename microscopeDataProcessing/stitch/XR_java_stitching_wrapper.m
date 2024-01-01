function [] = XR_java_stitching_wrapper(imageDirName, imageListFileName, varargin)
% wrapper for java stitching pipeline. 
% 
% 
% Inputs :   
%        imageDirName : full path for the directory of the images to be stitched
%   imageListFileName : full path for the coordinate information csv file
%
% Options (as 'specifier'-value pairs): 
%
%         'axisOrder' : Axis mapping of the coordinate system. Default: 'xyz'.
%         'ffcorrect' : flat field correction. Default: false.
%        'Resolution' : Image resolution in um, 2D (same for xy) or 3D vector.  
%         'resultDir' : name of stitching result directory. The directory is child directory of imageDirname
%      'parseCluster' : Use slurm-based cluster computing. Default: true. 
%         'jobLogDir' : Log directory for the slurm jobs. 
%       'cpusPerTask' : Number of cpus for a job. Default: 12
%         'Save16bit' : true|{false}. Save final results as 16bit or single. 
%              'uuid' : unique string for a job for saving files. 
%       'maxTrialNum' : Max number of times to rerun failure cases. 
%      'unitWaitTime' : For computing without cluster, the wait time per      
%                       file in minutes, in order to check whether the computing is done. 
%
%
% Author: Xiongtao Ruan (02/18/2020)
% 

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('imageDirName', @isstr);
ip.addRequired('imageListFileName', @isstr);
% ip.addParameter('Overwrite', true, @islogical);
ip.addParameter('axisOrder', 'x,z,-y', @isstr);
ip.addParameter('ffcorrect', false, @islogical);
ip.addParameter('Resolution', [], @isnumeric);
ip.addParameter('resultDir', 'stitching', @isstr);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('jobLogDir', '../job_logs', @isstr);
ip.addParameter('cpusPerTask', 12, @isnumeric);
ip.addParameter('Save16bit', false, @islogical);
ip.addParameter('uuid', '', @isstr);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 2, @isnumeric);

ip.parse(imageDirName, imageListFileName, varargin{:});

pr = ip.Results;
% Overwrite = pr.Overwrite;
axisOrder = pr.axisOrder;
ffcorrect = pr.ffcorrect;
Resolution = pr.Resolution;
resultDir = pr.resultDir;
jobLogDir = pr.jobLogDir;
parseCluster = pr.parseCluster;
cpusPerTask = pr.cpusPerTask;
Save16bit = pr.Save16bit;
uuid = pr.uuid;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;

% make root directory
stitching_rt = [imageDirName, filesep, resultDir];
mkdir(stitching_rt);

% temporary directory for intermediate results
stitching_tmp = [stitching_rt, filesep, 'tmp'];
mkdir(stitching_tmp);

% check if a slurm-based computing cluster exist
if parseCluster 
    [status, ~] = system('sinfo');
    if status ~= 0
        warning('A slurm-based computing cluster is not exist. Set parseCluster as false.')
        parseCluster = false;
    end
    if parseCluster && ~exist(jobLogDir, 'dir')
        warning('The job log directory does not exist, use ${stitching_tmp}/job_logs as job log directory')
        jobLogDir = sprintf('%s/job_logs', stitching_tmp);
        mkdir(jobLogDir);
    end
end

% uuid for the job
if isempty(uuid)
    uuid = get_uuid();
end

ffcorrect_str = string(ffcorrect);
Save16bit_str = string(Save16bit);

if numel(Resolution) == 2
    res_str = sprintf('%g,%g,%g', Resolution(1), Resolution(1), Resolution(2));
elseif numel(Resolution) == 3
    res_str = sprintf('%g,%g,%g', Resolution(1), Resolution(2), Resolution(3));    
end

%% parse image list information
% read image list csv file
t = readtable(imageListFileName, 'Delimiter','comma');
t_column_name = t.Properties.VariableNames;
fn = t.Filename;

specifyCam = true;
if all(~cellfun(@isempty, regexp(fn, '_Cam\w_ch', 'match')))
    expression = '(?<prefix>\w*)Scan_Iter_(?<Iter>\d+)(?<subIter>_?(\d+_)*\d+?)_Cam(?<Cam>\w+)_ch(?<ch>\d+)_CAM1_stack(?<stack>\d+)_(?<laser>\d+)nm_(?<abstime>\d+)msec_(?<fpgatime>\d+)msecAbs_(?<x>\d+)x_(?<y>\d+)y_(?<z>\d+)z_(?<t>\d+)t(?<suffix>_?\w*).tif';
elseif all(~cellfun(@isempty, regexp(fn, '_ch[0-9]_', 'match')))
    expression = '(?<prefix>\w*)Scan_Iter_(?<Iter>\d+)(?<subIter>_?(\d+_)*\d+?)_ch(?<ch>\d+)_CAM1_stack(?<stack>\d+)_(?<laser>\d+)nm_(?<abstime>\d+)msec_(?<fpgatime>\d+)msecAbs_(?<x>\d+)x_(?<y>\d+)y_(?<z>\d+)z_(?<t>\d+)t(?<suffix>_?\w*).tif';
    specifyCam = false;
end

tmp = regexpi(fn, expression, 'names');

matched_inds = true(numel(tmp), 1);

for f = 1:numel(tmp)
    if isempty(tmp{f})
        matched_inds(f) = false;
        continue;
    end
    
    t.prefix{f} = tmp{f}.prefix;
    t.Iter(f) = str2double(tmp{f}.Iter);
    t.subIter{f} = tmp{f}.subIter;
    t.fullIter{f} = [tmp{f}.Iter, tmp{f}.subIter];
    if specifyCam
        t.camera(f) = (tmp{f}.Cam);
    else
        % use N to represent cam if it is not contained in the filename
        t.camera(f) = 'N';
    end
    t.ch(f) = str2double(tmp{f}.ch);
    t.stack(f) = str2double(tmp{f}.stack);
    t.laser(f) = str2double(tmp{f}.laser);
    t.abstime(f) = str2double(tmp{f}.abstime);
    t.fpgatime(f) = str2double(tmp{f}.fpgatime);
    t.x(f) = str2double(tmp{f}.x);
    t.y(f) = str2double(tmp{f}.y);
    t.z(f) = str2double(tmp{f}.z);
    t.t(f) = str2double(tmp{f}.t);
end

t = t(matched_inds, :);

prefix = unique(t.prefix);
if ~isempty(prefix)
    prefix = prefix{1};
else
    prefix = '';
end
fullIter = unique(t.fullIter);
Ch = unique(t.ch);
Cam = unique(t.camera);
stackn = unique(t.stack);

% check whether the image files in the image list file exist 
dir_info = dir([imageDirName, filesep, '*.tif']);
imageFnames = {dir_info.name}';
image_file_exist_flag = true(numel(t.Filename), 1);
for f = 1 : numel(t.Filename)
    if ~contains(imageFnames, t.Filename{f})
        image_file_exist_flag(f) = false;
    end
end
if ~all(image_file_exist_flag)
    warning('Some files in the image list file do not exist! Ignore them in the stitching')
    disp(t.Filename(~image_file_exist_flag));
    t(~image_file_exist_flag, :) = [];
end


%% do stitching computing
row_exist_flag = true(numel(fullIter), numel(Cam), numel(Ch), numel(stackn)); % flag for whether the run exists.
is_done_flag = false(numel(fullIter), numel(Cam), numel(Ch), numel(stackn));
trial_counter = zeros(numel(fullIter), numel(Cam), numel(Ch), numel(stackn));
max_trial_num = maxTrialNum;

if parseCluster
    job_ids = -ones(numel(fullIter), numel(Cam), numel(Ch), numel(stackn));
    job_status_flag = false(numel(fullIter), numel(Cam), numel(Ch), numel(stackn));
end

while ~all(is_done_flag | trial_counter >= max_trial_num, 'all') ...
        || (parseCluster && any(job_status_flag & ~is_done_flag, 'all'))
    for n = 1:numel(fullIter)
        for ncam = 1:numel(Cam)
            for c = 1:numel(Ch)
                for s = 1:numel(stackn)
                    if is_done_flag(n, ncam, c, s) || trial_counter(n, ncam, c, s) >= max_trial_num
                        continue;
                    end
                    if ~row_exist_flag(n, ncam, c, s)
                        is_done_flag(n, ncam, c, s) = true;
                        continue;
                    end
                    task_id = sub2ind([numel(fullIter), numel(Cam), numel(Ch), numel(stackn)], n, ncam, c, s);

                    cur_t = t(t.ch == Ch(c) & t.camera == Cam(ncam) & strcmp(t.fullIter, fullIter{n}) & t.stack == stackn(s), :);

                    % obtain filenames                    
                    if isempty(cur_t)
                        row_exist_flag(n, ncam, c, s) = false;
                        continue;
                    end
                    
                    laser = unique(cur_t.laser);
                    abstime = unique(cur_t.abstime);
                    reltime = unique(cur_t.t);

                    % generate csv file for current run run
                    cur_ImageListDir = sprintf('%s/Iter_%s/Cam%s/ch%d/stack%04d', stitching_tmp, fullIter{n}, Cam(ncam), Ch(c), stackn(s));
                    mkdir_recursive(cur_ImageListDir);
                    cur_csv_fname = sprintf('%s/ImageList_Iter_%s_Cam%s_ch%d_stack%04d.csv', cur_ImageListDir, fullIter{n}, Cam(ncam), Ch(c), stackn(s));
                    writetable(cur_t(:, t_column_name), cur_csv_fname);

                    % also use flag based check of completion, to support
                    % distributed computing with same submission
                    stitch_save_dir = sprintf('%s/Cam%s/ch%d/', stitching_rt, Cam(ncam), Ch(c));
                    mkdir_recursive(stitch_save_dir);
                    stitch_save_fname = sprintf('%s/Scan_Iter_%s_Cam%s_ch%d_CAM1_stack%04d_%dnm_%06dmsec_%04dt.tif', ...
                        stitch_save_dir, fullIter{n}, Cam(ncam), Ch(c), stackn(s), laser, abstime, reltime);
                    cur_tmp_fname = sprintf('%s/Scan_Iter_%s_Cam%s_ch%d_CAM1_stack%04d_%dnm_%06dmsec_%04dt.tmp', ...
                        stitch_save_dir, fullIter{n}, Cam(ncam), Ch(c), stackn(s), laser, abstime, reltime);
                    
                    % first check if the computing is done
                    if exist(stitch_save_fname, 'file')
                        is_done_flag(n, ncam, c, s) = true; 
                        if exist(cur_tmp_fname, 'file')
                            delete(cur_tmp_fname);
                        end
                        continue;
                    end
                    
                    if exist(cur_tmp_fname, 'file')
                        % for cluster computing with master, check whether
                        % the job still alive. Otherwise, use waiting time
                        % for the check
                        if parseCluster
                            if job_ids(n, ncam, c, s) >= 0
                                job_status = check_slurm_job_status(job_ids(n, ncam, c, s), task_id);
                                % if job is still computing or pending
                                job_status_flag(n, ncam, c, s) = job_status >= 0;
                                % in case the old job is just finishing,
                                % also check whether the result already finished
                                if job_status >= 0 || exist(stitch_save_fname, 'file')
                                    continue;
                                end
                            end
                        else
                            per_file_wait_time = unitWaitTime; % minite
                            temp_file_info = dir(cur_tmp_fname);
                            if (datenum(clock) - [temp_file_info.datenum]) * 24 * 60 < numel(cur_t) * per_file_wait_time
                                continue; 
                            else
                                fclose(fopen(cur_tmp_fname, 'w'));
                            end
                        end
                    else
                        fclose(fopen(cur_tmp_fname, 'w'));
                    end

                    % stitching commands
                    if parseCluster
                        job_log_fname = [jobLogDir, '/job_%A_%a.out'];
                        job_log_error_fname = [jobLogDir, '/job_%A_%a.err'];
                                                
                        matlab_cmd = sprintf('setup([],true);java_stitching_frame_wrapper(''%s'',''%s'',''axisOrder'',''%s'',''ffcorrect'',%s,''Resolution'',[%s],''wavelength'',%d,''stitchResultFname'',''%s'',''Save16bit'',%s);toc;', ...
                            imageDirName, cur_csv_fname, axisOrder, ffcorrect_str, res_str, laser, stitch_save_fname, Save16bit_str);
                        stitch_cmd = sprintf('module load matlab/r2021a; matlab -nodisplay -nosplash -nodesktop -r \\"%s\\"', matlab_cmd);
                        cmd = sprintf('sbatch --array=%d -o %s -e %s -p abc --qos abc_normal -n1 --mem-per-cpu=20G --cpus-per-task=%d --wrap="%s"', ...
                            task_id, job_log_fname, job_log_error_fname, cpusPerTask, stitch_cmd);
                        [status, cmdout] = system(cmd, '-echo');

                        job_id = regexp(cmdout, 'Submitted batch job (\d+)\n', 'tokens');
                        job_id = str2double(job_id{1}{1});
                        job_ids(n, ncam, c, s) = job_id;
                    else
                        java_stitching_frame_wrapper(imageDirName, cur_csv_fname, 'axisOrder', axisOrder, ...
                            'ffcorrect', ffcorrect, 'Resolution', Resolution, 'wavelength', laser, ...
                            'stitchResultFname', stitch_save_fname, 'Save16bit', Save16bit, 'uuid', uuid);
                    end
                    
                    trial_counter(n, ncam, c, s) = trial_counter(n, ncam, c, s) + 1;    
                    
                    % check if computing is done
                    if exist(stitch_save_fname, 'file')
                        is_done_flag(n, ncam, c, s) = true;
                        if exist(cur_tmp_fname, 'file')
                            delete(cur_tmp_fname);
                        end
                    end
                end
            end
        end
    end
    
    % wait 30 seconds if some tasks are still computing
    if ~all(is_done_flag | trial_counter >= max_trial_num, 'all') || ...
            (parseCluster && any(job_status_flag & ~is_done_flag, 'all'))
        pause(30);
    end
end

end

