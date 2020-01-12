function mismatchRgbMergeBatchSetup(inputDir, varargin)
%This program reads through a folder containing two exposures of
%vimentin in folders with 'Dim' and 'Bright' in their title and loop
%through the folder to create hdrMerge images in their own folder for each
%one.
%-Jessica Tytell December 5, 2011

% parse input
ip = inputParser;
ip.addRequired('inputDir', @ischar);
ip.addOptional('primeChan', 'images', @ischar);
ip.addOptional('secChan', ['hdrMerge' filesep 'logMysteryTIFFs'], @ischar); 
ip.parse(inputDir, varargin{:});
primeChan = ip.Results.primeChan;
secChan = ip.Results.secChan;

%Get the folders for each movie
movieFolders = dir([inputDir filesep '*EB*VIM*']);
movieFolders = movieFolders(arrayfun(@(x)(x.isdir && ... %Retain only the directories. Do it this way so it works on linux and PC
    ~(strcmp(x.name,'.') || strcmp(x.name,'..'))),movieFolders));
 

nMovies = length(movieFolders);

%loop through movies
for j = 1:nMovies
    
    disp(['Processing folder ' num2str(j) ' of ' num2str(nMovies)])
    
    %Get current folder path for readability
    currDir = [inputDir filesep movieFolders(j).name];
    disp(currDir);
    
    %exit loop if files already exist - no overwrite allowed.
    if exist([currDir filesep 'VimentinPlusTipMerge'], 'dir')
        disp('Vimentin Merge file already exists. No override function exists, please delete file and start again');
        continue;
    else
        %get base and intermittent file directories
        baseDir = [currDir filesep primeChan];
        intermitDir = [currDir filesep secChan];
        disp(baseDir);
        disp(intermitDir);
        
        %exit loop with warning if no images in directory
        if isempty(baseDir)
            disp('Main image folder is empty: Please try again');
        elseif isempty(intermitDir)
            disp('Secondary image folder is empty: Please try again');
        else
            
            %send to hdrMergeFileBatch
            mismatchRgbMergeBatch(baseDir, intermitDir);
        end
        
    end
end
