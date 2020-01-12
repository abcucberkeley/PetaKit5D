function hdrMergeFileBatchSetup(inputDir)
%This program reads through a folder containing two exposures of
%vimentin in folders with 'Dim' and 'Bright' in their title and loop
%through the folder to create hdrMerge images in their own folder for each
%one.
%-Jessica Tytell December 5, 2011


% parse input
ip = inputParser;
ip.addRequired('inputDir', @ischar);
ip.parse(inputDir);

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
    if exist([currDir filesep 'hdrMerge'], 'dir')
        disp('hdr file already exists. No override function exists, please delete file and start again');
        continue;
    else
        %get bright and dim file directories
        dimStruct = dir([currDir filesep '*Dim*']);
        brightStruct = dir([currDir filesep '*Bright*']);
        disp(dimStruct);
        disp(brightStruct);
        
        %get directory names to pass to batch function
        if length(dimStruct) > 1
            disp('Too many folders labeled "Dim" ');
        else
            dimName = dimStruct.name;
            dimDir = [currDir filesep dimName];
        end
        
        if length(brightStruct) > 1
            disp('Too many folders labeled "Bright"');
        else
            brightName = brightStruct.name;
            brightDir = [currDir filesep brightName];
        end
        
        %send to hdrMergeFileBatch
        hdrMergeFileBatch(dimDir, brightDir);
        
    end
    
end
