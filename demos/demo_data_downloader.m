function [] = demo_data_downloader(destPath, dataType)
% automatically download demo dataset from Dropbox
% dataType: light_sheet or other_modalities


if nargin < 2
    dataType = 'light_sheet';
end

switch lower(dataType)
    case 'light_sheet'
        fprintf('Download cell data with light sheet microscopy...\n')
        demo_light_sheet_cell_data_downloader(destPath);
    case 'other_modalities'
        fprintf('Download data from other imaging modelities (2-photon, Confocal, Phase, and Widefield)...\n')
        demo_other_modalities_data_downloader(destPath);
    otherwise 
        error('Unknown data type. Data type must be light_sheet or other_modelities!')
end


end

