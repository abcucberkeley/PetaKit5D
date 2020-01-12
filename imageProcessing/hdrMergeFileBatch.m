function hdrMergeFileBatch(dimDir, brightDir)
%This function takes an input of directories containing the bright and dim images
%needed for high dynamic range merging of 16bit images. It then calls the
%highDynamicRangeMerge function on each and writes the output to nested folders inside a new
%folder called hdrMergedImages.

%-Jessica Tytell December 5, 2011

% parse input
ip = inputParser;
ip.addRequired('dimDir', @ischar); 
ip.addRequired('brightDir', @ischar);


%read in files to arrays
dimFiles = imDir(dimDir);
brightFiles = imDir(brightDir);

%set location and make output dir
outputLocation = fileparts(dimDir);
outputDir = ([outputLocation filesep 'hdrMerge']);
mkdir(outputDir);

%show files to test reading - comment out for final version
% disp('dimFiles = ' );
% disp(dimFiles);
% disp('brightFiles = ' );
% disp(brightFiles);

dimNames = {dimFiles.name};
brightNames = {brightFiles.name};

%check to make sure there are the same # of files
nDim = length(dimNames);
nBright = length(brightNames);
if nDim == nBright
    disp(['There are ', num2str(nDim), ' files in both folders. Win!']);
    %make new directory for output images
    %go up one directory
%     upDir = [dimDir filesep '..'];
%     cd(upDir);
%     make new directory for merged images
%    mkdir('hdrMergedImages');
%     cd('hdrMergedImages');
    mkdir([outputDir filesep 'combinedMysteryScaledTIFFs']);
    mkdir([outputDir filesep 'logMysteryTIFFs']);
    mkdir([outputDir filesep 'logImageScaledTIFFs']);
    
    %establish waitbar
    h = waitbar(0,'<jeopardy music>');
    
    %loop through images    
    for j = 1:nDim
        nextDimFile = [dimDir filesep dimNames{j}];
        nextBrightFile = [brightDir filesep brightNames{j}];
        [~, ~, time] = getFilenameBody(dimNames{j});
        dimIm = double(imread(nextDimFile));
        brightIm = double(imread(nextBrightFile));
        [combinedMystery, logMystery, logImage] = highDynamicRangeMerge(dimIm, brightIm);
        
        %update waitbar
         waitbar(j / nDim)
        %display images that are returned
%         figure('Name', 'CombinedMystery'), imshow(combinedMystery, []);
%         figure('Name', 'logMystery'), imshow(logMystery, []);
%         figure('Name', 'logImage'), imshow(logImage, []);
%         
        %attempt to write images
        maxNumMyst = max(combinedMystery(:));
        combinedMysteryScaled = scaleContrast(combinedMystery, [0 maxNumMyst], [0 255]);
        outputName = [outputDir filesep 'combinedMysteryScaledTIFFs' filesep 'hdrMerge_combMyst_' time '.TIF'];
        imwrite(uint8(combinedMysteryScaled), outputName, 'tiff');
        
        maxNumLog = max(logMystery(:));
        minNumLog = min(logMystery(:));
        logMysteryScaled = scaleContrast(logMystery, [minNumLog maxNumLog], [0 255]);
        outputNameLog= [outputDir filesep 'logMysteryTIFFs' filesep 'hdrMerge_logMyst_', time '.TIF'];
        imwrite(uint8(logMysteryScaled), outputNameLog, 'tiff');
        
        maxNumLogIm = max(logImage(:));
        minNumLogIm = min(logImage(:));
        logImageScaled = scaleContrast(logImage, [minNumLogIm maxNumLogIm], [0 255]);
        outputNameLogIm= [outputDir filesep 'logImageScaledTIFFs' filesep 'hdrMerge_logImage_', time '.TIF'];
        imwrite(uint8(logImageScaled), outputNameLogIm, 'tiff');
%         imtool(logImageScaled,[]);
    end
    close (h);
    
else disp(['There are ', nDim, ' files in the Dim channel and ', nBright, ' images in the bright folder. Image folders must contain the same number of images. Input FAIL.']);
end
