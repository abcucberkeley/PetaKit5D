function mismatchRgbMergeBatch(baseFileDir, intermitDir)
%This function is designed to make RGB images of two channels where there
%are different numbers of time points for each. It now uses Sylvain's
%KDtree function to find the closest secondary channel to assign the
%primary channel to and then loops through using these indeces to make
%rgbMovie

%-Jessica Tytell December 2, 2011 
% Sebastien Besson (Assist) December 6, 2011

% parse input
ip = inputParser;
ip.addRequired('baseFileDir', @ischar); 
ip.addRequired('intermitDir', @ischar);
ip.parse(baseFileDir, intermitDir);

%set location and make output dir
outputLocation = fileparts(baseFileDir);
outputDir = ([outputLocation filesep 'VimentinPlusTipMerge']);
mkdir(outputDir);

%read in files and show in window
[baseFiles,~,baseTimes] = imDir(baseFileDir);
% disp('baseFiles = ' );
% disp(baseFiles);
[intermitFiles,~,intermitTimes] = imDir(intermitDir);
% disp('intermittent Files = ');
% disp(intermitFiles);

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
    
    
    
     