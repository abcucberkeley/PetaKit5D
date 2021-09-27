function [] = separatePhases(dataFile, varargin)

%fn = '/clusterfs/fiona/Data/20210923_latticeSIM/data06_100perc/RAW_488_slow_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0013199486msecAbs_000x_000y_000z_0000t.tif';
%dataFile = '/clusterfs/fiona/Data/20210923_latticeSIM/data06_100perc/RAW_exp08_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0006674367msecAbs_000x_000y_000z_0000t.tif';
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataFile');
ip.addParameter('nphases', 5, @isnumeric);

ip.parse(dataFile, varargin{:});

pr = ip.Results;
nphases = pr.nphases;

data = readtiff(dataFile);
for p=1:nphases
    writetiff(data(:,:,p:nphases:end),[dataFile(1:end-4) '_phase' num2str(p) '.tif']);
end

end
