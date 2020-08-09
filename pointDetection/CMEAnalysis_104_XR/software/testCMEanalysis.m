% Now using MATLAB-based testing/performance suite
% Test cmeAnalysis
% Using  locally stored data for now and only testing analysis portion, not detectino/tracking
% Andrew R. Jamieson 2016
%
% Copyright (C) 2017, Danuser Lab - UTSouthwestern 
%
% This file is part of CMEAnalysis_Package.
% 
% CMEAnalysis_Package is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% CMEAnalysis_Package is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with CMEAnalysis_Package.  If not, see <http://www.gnu.org/licenses/>.
% 
% 


% function scriptCMEautotest()

% preconditions
disp(['PWD:', pwd]);
% Dump path for debugging
s_path = strsplit(path,':');
s_match = (cellfun(@(x) regexp(x,'toolbox'), s_path, 'UniformOutput', false))';
matlab_paths = s_path(cellfun(@isempty, s_match))';
disp('    [MATLAB] current top-level paths....');
disp(matlab_paths);

disp(['Java heap max: ' num2str(java.lang.Runtime.getRuntime.maxMemory/1e9) 'GB'])
disp('Starting cmeAnalysis test script');

% ----------------------------------------------------------------------------
skip_detection_N_tracking = false
small_data_set = true

if skip_detection_N_tracking
	% Quick testing--- without running detection and tracking....
	if small_data_set
        disp('Running on truncated data set - small');
		if strcmp(computer('arch'),'win64')
			data_root = 'C:\Users\Andrew\Data\raw\CME\xinxin\auto_test\small\';
			data_rootEpi = 'C:\Users\Andrew\Data\raw\CME\xinxin\auto_test\';
		else
			data_root = '/work/bioinformatics/s170480/Data/CME/xinxin/auto_test/small/';
			data_rootEpi = '/work/bioinformatics/s170480/Data/CME/xinxin/auto_test/';
		end			
    else
        disp('Running on large data set - may take some time...');
        if strcmp(computer('arch'),'win64')
			data_root = 'C:\Users\Andrew\Data\raw\CME\xinxin\auto_test\';
		else
			data_root = '/work/bioinformatics/s170480/Data/CME/xinxin/auto_test/';
		end
		data_rootEpi = 	data_root;
	end

	% cmeData_dir = fullfile(data_root, ['single_channel' filesep 'Zuzana_ARPE_CLC_egfp_07032014_control']);
	cmeData_dir2 = fullfile(data_root, ['single_channel' filesep 'Zuzana_ARPE_CLC_egfp_07032014_mu_PIP2_mutant']);
	% epiTIRFData_dir = fullfile(data_root, ['epi_tirf' filesep 'Zuzana_arpe_clc_egfp_14082015_control']);
	epiTIRFData_dir2 = fullfile(data_rootEpi, ['epi_tirf' filesep 'Zuzana_arpe_clc_egfp_14082015_mu_PIP2_mutant']);
	% ----------------------------------------------------------------------------
else
	% - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% ----------------------------------------------------------------------------
	disp('Testing on "clean" data -- including detection and tracking')
	% no mask yet.. and no tracking
	if small_data_set
		disp('Running on truncated data set - small');
        if strcmp(computer('arch'),'win64')
			data_root = 'C:\Users\Andrew\Data\raw\CME\xinxin\auto_test\clean\small\';
			data_rootEpi = 'C:\Users\Andrew\Data\raw\CME\xinxin\auto_test\clean\epi_tirf\';
		else
			data_root = '/work/bioinformatics/s170480/Data/CME/xinxin/auto_test/clean/small/';
			data_rootEpi = '/work/bioinformatics/s170480/Data/CME/xinxin/auto_test/clean/epi_tirf/';
		end			
    else
        disp('Running on large data set - may take some time...');
        if strcmp(computer('arch'),'win64')
			data_root = 'C:\Users\Andrew\Data\raw\CME\xinxin\auto_test\clean\';
		else
			data_root = '/work/bioinformatics/s170480/Data/CME/xinxin/auto_test/clean/';
		end
		data_rootEpi = 	data_root;
	end

	% ----Initialization of temp dir
	package_name = 'cmeAnalysis';
	t_stamp = datestr(now,'ddmmmyyyyHHMMSS');
	tmpdir = fullfile(tempdir, [package_name '_test_' t_stamp]);
	mkdir(tmpdir);

	tic; 
	disp(['Copying ' data_root 'to ' tmpdir]);
	copyfile(data_root, tmpdir);
    if ~strcmp(data_rootEpi,data_root)
        disp('copying epi tirf data also...')
        copyfile(data_rootEpi, [tmpdir filesep 'epi_tirf' filesep]);
    end
	disp('Done copying')
	toc

	% cmeData_dir = fullfile(data_root, ['single_channel' filesep 'Zuzana_ARPE_CLC_egfp_07032014_control']);
	cmeData_dir2 = fullfile(tmpdir, ['single_channel' filesep 'Zuzana_ARPE_CLC_egfp_07032014_mu_PIP2_mutant']);

	% epiTIRFData_dir = fullfile(data_root, ['epi_tirf' filesep 'Zuzana_arpe_clc_egfp_14082015_control']);
	epiTIRFData_dir2 = fullfile(tmpdir, ['epi_tirf' filesep 'Zuzana_arpe_clc_egfp_14082015_mu_PIP2_mutant']);
	% ----------------------------------------------------------------------------
end
% - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sortTiffStacks; % helper function for organizing data structure

% Setup data to use loadConditionData
% data = loadConditionData(cmeData_dir, {''}, {'EGFP'}); %[works]
data2 = loadConditionData(cmeData_dir2, {''}, {'EGFP'}); %[works]

%% RUN__cmeAnalysis__Data1
% results = cmeAnalysis(data);  %[works]
results2 = cmeAnalysis(data2);  %[works]
plotCellArea({results2.lftRes}, {'ctrl'}); % (included in cmeAnalysis)
plotInitiationDensity({results2.lftRes, results2.lftRes}, {'ctrl', 'other'}) % error bar styling MATLAB version update

% ---
%  cmeAnalysis(data, 'Overwrite', [true false false true]) -- manually go
%  over tracks for semgnetation mask.
% ---

%% Run ccpSorter
% ccpSorter(data);
ccpSorter(data2);

%% Run cmeDataViewer
% cmeDataViewer(data(1));
% cmeDataViewer(data(2)); 
cmeDataViewer(data2(1));
cmeDataViewer(data2(2));

%% Run analyzeBleaching
% analyzeBleaching(data(1));
% analyzeBleaching(data(2));
analyzeBleaching(data2(1));
analyzeBleaching(data2(2));

%% RUN__runLifetimeAnalysis__and__Plots
% [lftRes, res] = runLifetimeAnalysis(data);  % error bar fixes
[lftRes2, res2] = runLifetimeAnalysis(data2);  % error bar fixes
% plotLifetimes(lftRes);  % works (included in cmeAnalysis)
plotLifetimes(lftRes2);  % works (included in cmeAnalysis)
plotLifetimeComparison({lftRes2, lftRes2}, {'con1','con2'}); % wrap with cell % ask Marcel? (need at least two data)

%% loadTracks__and__plots
track_num = 1;
movie_num = 1;
track = loadTracks(data2(movie_num)); % xinxin 
plotTrack(data2(movie_num), track(track_num)); %  % (included in cmeAnalysis)
plotIntensityCohorts(data2); % (included in cmeAnalysis)
plotMaxIntensityDistribution(data2);  % works, but old % MATLAB update fixes
plotIntensityDistributions(data2(movie_num)); % works
% ---
% plotInitIntensityVsLifetime(data);  % ???? not used regularly ????
% ---

%% epiTIRFAnalysis (with prev tracking)
dataEpi2 = loadConditionData(epiTIRFData_dir2, {'TIRF', 'EPI'}, {'EGFP', 'EGFP'});
epiTIRFAnalysis({dataEpi2, dataEpi2},  {'EGFP', 'EGFP'}, 9.5); % Still in development
