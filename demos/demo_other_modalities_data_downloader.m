function [] = demo_other_modalities_data_downloader(destPath)
% automatically download demo dataset from Dropbox


if nargin == 0 || isempty(destPath)
    warning('The destination path does not exist, save dataset to ~/Downloads/!')
    if ispc
        destPath = fullfile(getenv('USERPROFILE'), 'Downloads');   
    else
        destPath = '~/Downloads/';
    end
end

if ispc
    destPath = strrep(destPath, '\', '/');
end

dataPath = [strip(destPath, 'right', '/'), '/PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets/'];

data_types = {'2P', 'Confocal', 'Phase', 'Widefield'};

% check if dataset already downloaded
if exist(dataPath, 'dir')
    % check if all necessary files exist
    dir_info = dir([dataPath, '/*']);
    fsns = {dir_info.name};
    
    file_exist = true;
    
    % readme.txt
    file_exist = file_exist && sum(matches(fsns, "readme.txt")) == 1;
    
    file_exist = file_exist && all(ismember(data_types, fsns));

    % check if the individual datasets exist
    if file_exist
        for i = 1 : numel(data_types)
            data_type = data_types{i};
            dataPath_i = sprintf('%s/%s/', dataPath, data_type);
            file_exist = check_single_dataset(dataPath_i, data_type);
            if ~file_exist
                break;
            end
        end
    end
    
    if file_exist
        fprintf('The datasets already exist in "%s", skip downloading!\n', dataPath);
        return;
    end
end

% The demo dataset is available in zenodo (https://zenodo.org/records/11500863). We also shared the dataset from
% Drobox to allow much faster downloads. 
% url = 'https://zenodo.org/records/11500863/files/PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets.tar?download=1';
url = 'https://www.dropbox.com/scl/fi/sp1pz423jw6gpozab5yrf/PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets.tar?rlkey=i4wrwvphm9hz78z3kcqomqm1b&st=az3w6z0v&dl=1';

fprintf('Download other modalities (2-photon, Confocal, Phase, and Widefield) demo datasets from Dropbox...\n')
filename = [destPath, '/PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets.tar'];
outputfilename = websave(filename, url);

fprintf('Download finished! \nUntar dataset...\n')
if ~ismac
    untar(outputfilename, destPath);
else
    system(sprintf('tar -vxf %s/PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets.tar -C %s', destPath, destPath));
end

fprintf('Untar dataset finished! \nDelete the tar file...\n');
delete(outputfilename);

fprintf('Done!\n');

end


function [file_exist] = check_single_dataset(dataPath, data_type)

file_exist = true;
switch data_type
    case '2P'
        dir_info = dir([dataPath, '/*']);
        fsns = {dir_info.name};

        % tile files
        file_exist = file_exist && sum(matches(fsns, "Tile" + wildcardPattern + "tif")) == 8;

        % Flat field
        file_exist = file_exist && sum(matches(fsns, "flatfield")) == 1;

        % check if flatfield files exist
        file_exist = file_exist && exist([dataPath, 'flatfield/background.tif'], 'file') && exist([dataPath, 'flatfield/flatfield.tif'], 'file');
    case 'Confocal'
        file_exist = file_exist && exist([dataPath, 'Confocal_PSF.tif'], 'file') && exist([dataPath, 'Confocal_raw.tif'], 'file');
    case 'Phase'
        dir_info = dir([dataPath, '/*']);
        fsns = {dir_info.name};

        % tile files
        file_exist = file_exist && sum(matches(fsns, "Scan" + wildcardPattern + "tif")) == 40;
        
        % image List file
        file_exist = file_exist && sum(matches(fsns, "ImageList_all_timepoints.csv")) == 1;

        % ff_corrected
        file_exist = file_exist && sum(matches(fsns, "ff_corrected")) == 1;

        if file_exist
            dir_info = dir([dataPath, 'ff_corrected/*tif']);
            fsns = {dir_info.name};
            file_exist = file_exist && sum(matches(fsns, "Scan" + wildcardPattern + "tif")) == 40;
        end
    case 'Widefield'
        file_exist = file_exist && exist([dataPath, 'WF_PSF.tif'], 'file') && exist([dataPath, 'WF_raw.tif'], 'file');
    otherwise 
        error('Unknown data type!');
end

end

