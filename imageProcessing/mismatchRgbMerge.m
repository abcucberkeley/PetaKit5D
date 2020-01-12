function mismatchRgbMerge()
%This function is designed to make RGB images of two channels where there
%are different numbers of time points for each. It reads in two directories
%and overlays the channel with less time points with the closest time point
%from the base frame using Sylvain's KDtree finder.

%-Jessica Tytell December 2, 2011
% Sebastien Besson (Assist) December 6, 2011 - adding new algorithm

%find directories containing dim and bright images
baseFileDir = uigetdir('', 'Select Folder containing images in each time point (base file)');
intermitDir = uigetdir('', 'Select Folder with intermittent images');
outputDir = uigetdir('', 'Select Output Folder for merged images');

%read in files and show in window
[baseFiles,~,baseTimes] = imDir(baseFileDir);
disp('baseFiles = ' );
disp(baseFiles);
[intermitFiles,~,intermitTimes] = imDir(intermitDir);
disp('intermittent Files = ');
disp(intermitFiles);

if isempty(baseTimes) || isempty(intermitTimes)
    error('One of the image directories is empty. Unfortunately this program does not work with theoretical data');
end

%write names to own array as string
baseNames = {baseFiles.name};
intermitNames = {intermitFiles.name};

%get length of longer file
nBase = length(baseNames);

%Find closest time index in the intermittent timepoints
intermitIndex=KDTreeClosestPoint(intermitTimes,baseTimes);

%add waitbar for impatient people (or people who don't trust this code)
h = waitbar(0,'Please wait...');

%set up incrementer for intermittent files
% inc = 1;

% Create formatted string for padding zeros
fString = ['%0' num2str(floor(log10(nBase))+1) '.f'];

for j = 1:nBase
    %get next file and read in
    nextBaseFile = [baseFileDir filesep baseNames{j}];
    baseIm = double(imread(nextBaseFile));
    
    % get time point from filename
    [~, baseBody] = getFilenameBody(baseNames{j});
    
    %test if files extensions are the same or smaller
%     [~, ~, intTime] = getFilenameBody(intermitNames{inc});
    
    %test to see if next intermittent time point is higher than current
    %base time
%     if str2double(baseTime) > str2double(intTime)
%         
%         inc = inc+1;
%         
%     end
     
    nextIntFile = [intermitDir filesep intermitNames{intermitIndex(j)}];
    intermitIm = double(imread(nextIntFile));
    
    %rgb merge
    rgbIm = ch2rgb(baseIm, intermitIm, []);
    
    outputName= [outputDir filesep baseBody '_merge' num2str(baseTimes(j),fString) '.TIF'];
    imwrite((rgbIm), outputName, 'tiff');
    
    %update waitbar
    waitbar(j / nBase)

end
close(h);

disp('<snoopy dance>');
    
    
    

% %read in files and show in window
% baseFiles = imDir(baseFileDir);
% disp('baseFiles = ' );
% disp(baseFiles);
% intermitFiles = imDir(intermitDir);
% disp('intermittent Files = ');
% disp(intermitFiles);
% 
% %write names to own array as string
% baseNames = {baseFiles.name};
% intermitNames = {intermitFiles.name};
% 
% %get length of longer file
% nBase = length(baseNames);
% 
% %add waitbar for impatient people (or people who don't trust this code)
% h = waitbar(0,'Please wait...');
% 
% %set up incrementer for intermittent files
% inc = 1;
% for j = 1:nBase
%     %get next file and read in
%     nextBaseFile = [baseFileDir filesep baseNames{j}];
%     baseIm = double(imread(nextBaseFile));
%     
%     % get time point from filename
%     [~, baseBody, baseTime] = getFilenameBody(baseNames{j});
%     
%     %test if files extensions are the same or smaller
%     [~, ~, intTime] = getFilenameBody(intermitNames{inc});
%     
%     %test to see if next intermittent time point is higher than current
%     %base time
%     if str2double(baseTime) > str2double(intTime)
%         
%         inc = inc+1;
%         
%     end
%     
%     nextIntFile = [intermitDir filesep intermitNames{inc}];
%     intermitIm = double(imread(nextIntFile));
%     
%     %rgb merge
%     rgbIm = ch2rgb(baseIm, intermitIm, []);
%     
%     outputName= [outputDir filesep baseBody '_merge' baseTime '.TIF'];
%     imwrite((rgbIm), outputName, 'tiff');
%     
%     %update waitbar
%     waitbar(j / nBase)
% 
% end
% close(h);
% 
% disp('<snoopy dance>');
    
    
    
    