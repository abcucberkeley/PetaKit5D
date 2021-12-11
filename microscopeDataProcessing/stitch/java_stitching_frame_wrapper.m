function java_stitching_frame_wrapper(imageDirName, imageListFileName, varargin)
% java stitching wrapper for each frame/channel. 
% 
% Author: Xiongtao Ruan (02/21/2020)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('imageDirName', @isstr);
ip.addRequired('imageListFileName', @isstr);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('axisOrder', 'x,z,-y', @isstr);
ip.addParameter('ffcorrect', false, @islogical);
ip.addParameter('Resolution', [], @isnumeric);
ip.addParameter('wavelength', [], @isnumeric);
ip.addParameter('stitchResultFname', '', @isstr);
ip.addParameter('Save16bit', false, @islogical);
ip.addParameter('uuid', '', @isstr);

ip.parse(imageDirName, imageListFileName, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
axisOrder = pr.axisOrder;
ffcorrect = pr.ffcorrect;
Resolution = pr.Resolution;
stitchResultFname = pr.stitchResultFname;
Save16bit = pr.Save16bit;
uuid = pr.uuid;
wavelength = pr.wavelength;

% check whether the result already exist
if exist(stitchResultFname, 'file') && ~Overwrite
    fprintf('Stitching result %s already exist! Skip the computing!\n', stitchResultFname);
    return;
end

if isempty(uuid)
    uuid = get_uuid();
end

if numel(Resolution) == 2
    res_str = sprintf('%g,%g,%g', Resolution(1), Resolution(1), Resolution(2));
elseif numel(Resolution) == 3
    res_str = sprintf('%g,%g,%g', Resolution(1), Resolution(2), Resolution(3));    
end

% add the parent folder to path.
% fpath = fileparts(which(mfilename));
% cd(fpath);
% addpath(genpath([fpath, filesep, '..']));

%% stitching
stitch_cmd_prefix = 'bash ./XiongtaoScripts/stitch/java_stitching_frame_commands.sh';
if ffcorrect
    stitch_cmd_prefix = [stitching_cmd_prefix, ' -f'];
end

stitch_cmd = sprintf('%s -a %s -b %s -c %d -i %s -r %s', stitch_cmd_prefix, ...
    axisOrder, imageDirName, wavelength, imageListFileName, res_str);

cmd = stitch_cmd;
[status, cmdout] = system(cmd, '-echo');


%% save stitched results as a single tiff file

fprintf('Collect stitching results...\n');

stitch_save_fname = stitchResultFname;

stitch_save_dir = fileparts(stitch_save_fname);
mkdir_recursive(stitch_save_dir);

ImageListDir = fileparts(imageListFileName);
stitch_result_dir = sprintf('%s/slice-tiff-s0/ch0/', ImageListDir);
dir_info = dir([stitch_result_dir, '*tif']);
if isempty(dir_info)
    return;
end

slice_fnames = {dir_info.name};

% the names are just numbers starting from 0
% first check tomake sure the files are complete
nz = numel(slice_fnames);

for z = 1 : nz
    fname = sprintf('%d.tif', z - 1);
    if ~any(strcmp(slice_fnames, fname))
        error('%s does not exist! Please check to make sure the stitching results are complete.', fname);
    end
end

img_z = readtiff([stitch_result_dir, slice_fnames{1}]);
img_stitch = zeros([size(img_z), nz]);
for z = 1 : nz
    fname = sprintf('%d.tif', z - 1);
    img_z = readtiff([stitch_result_dir, fname]);
    img_stitch(:, :, z) = img_z;
end

% avoid different jobs save the same file
stitch_save_fname_tmp = sprintf('%s%s%s', stitch_save_fname(1 : end - 4), uuid, stitch_save_fname(end - 3 : end));
if Save16bit
    writetiff(uint16(img_stitch), stitch_save_fname_tmp);
else
    writetiff(single(img_stitch), stitch_save_fname_tmp);
end
movefile(stitch_save_fname_tmp, stitch_save_fname);

end

