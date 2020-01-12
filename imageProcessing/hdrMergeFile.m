function hdrMergeFile()
%This function asks for directories containing the bright and dim images
%needed for high dynamic range merging of 16bit images. It then calls the
%highDynamicRangeMerge function on each and writes the output to a new
%folder called hdrMergedImages.
%Note, down the road, I would like to change this program to be able to
%change the output directory as well and to include the options to change
%the mystery factor # (assuming that it hasn't been canceled by then).
%-Jessica Tytell October 25, 2011

%find directories containing dim and bright images
dimDir = uigetdir('', 'Select Folder containing Dim images');
brightDir = uigetdir('', 'Select Folder containing Bright images');
outputDir = uigetdir('', 'Select Output Folder');
% %for troubleshooting only
% dimDir = '/home/jdt2/Desktop/IFproject/110915/D1_006/w2491_VimDim';
% brightDir = '/home/jdt2/Desktop/IFproject/110915/D1_006/w3491_VimBright';

%find parent directory containing both images
% THIS SECTION IS NOT WORKING - something wrong wigh regexp matching....
%         inputDir = uigetdir('', 'Select Folder containing dim and bright images in separate folders');
%         dimDir = dirList(cellfun(@any,regexp({inputDir},'Dim')));
%         brightDir = dirList(cellfun(@any,regexp({inputDir},'Bright')));


%read in files to arrays
dimFiles = imDir(dimDir);
brightFiles = imDir(brightDir);

disp('dimFiles = ' );
disp(dimFiles);
disp('brightFiles = ' );
disp(brightFiles);

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

