function movieData = backgroundSubtractMovie(movieData,paramsIn)
%BACKGROUNDSUBTRACTMOVIE corrects for background by subtracting the average background value
%
% movieData = backgroundSubtractMovie(movieData)
%
% movieData = backgroundSubtractMovie(movieData,paramsIn)
%
% This function performs background subtraction on the movie described by
% the input movieData. This is accomplished by averaging the intensity in
% the areas covered by background masks, as created using
% createMovieBackgroundMasks.m. This average value is determined for each
% frame and then subtracted from each pixel in that frame. Negative values
% are converted to zero. 
% 
% 
% Input:
% 
%   movieData - The MovieData object describing the movie, as created using
%   setupMovieDataGUI.m
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below:
% 
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('OutputDirectory' -> character string) Optional. A character
%       string specifying the directory to save the corrected images to.
%       Corrected images for different channels will be saved as
%       sub-directories of this directory. If not input, the corrected
%       images will be saved to the same directory as the movieData, in a
%       sub-directory called "background_subtracted_images"
%
%       ('ChannelIndex'-> Positive integer scalar or vector) The integer
%       indices of the channel(s) to perform shade correction on. This
%       index corresponds to the channel's location in the array
%       movieData.channels_. If not input, all channels will be corrected.
%      
%       ('MaskChannelIndex' -> Positive integer scalar or vector)
%       This parameter specifies the channels to use masks from when
%       performing background subtraction on each channel. This allows
%       masks generated for one channel to be used in performing background
%       subtraction on other channels, which may or may not themselves have
%       masks. This vector or scalar must be the same size as ChannelIndex.
%       Optional. If not input, each channel will use it's own masks. 
%
%       ('BatchMode' -> True/False)
%       If this option value is set to true, all graphical output and user
%       interaction is suppressed.
%
%
%
% Output:
%
%   movieData - the updated movieData object with the correction
%   parameters, paths etc. stored in it, in the field movieData.processes_.
%
%   The corrected images are written to the directory specified by the
%   parameter OuptuDirectory, with each channel in a separate
%   sub-directory. They will be stored as bit-packed .tif files. 
%
% 
% Hunter Elliott, 11/2009
% Revamped 5/2010
%

%%  --------- Parameters ------- %%

pString = 'bs_'; %The string to prepend before the background-subtracted image directory & channel name
saveName = 'background_subtraction_values_for_channel_'; %File name for saving subtracted values
dName = 'background_subtracted_images_for_channel_';%String for naming the directories for each corrected channel

%% ----------- Input ------------ %%



%Check that input object is a valid moviedata
if ~isa(movieData,'MovieData')
    error('The first input argument must be a valid MovieData object!')
end

if nargin < 2
    paramsIn = [];
end

%Get the indices of any previous background subtraction processes from this
%function
iProc = find(cellfun(@(x)(isa(x,'BackgroundSubtractionProcess')),movieData.processes_),1);                          

%If the process doesn't exist, create it with default settings.
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(BackgroundSubtractionProcess(movieData,movieData.outputDirectory_));                                                                                                 
end

%Parse input, store in parameter structure
p = parseProcessParams(movieData.processes_{iProc},paramsIn);

nChan = numel(movieData.channels_);


if max(p.ChannelIndex) > nChan || min(p.ChannelIndex)<1 ...
        || ~isequal(round(p.ChannelIndex),p.ChannelIndex)
    error('Invalid channel numbers specified! Check ChannelIndex input!!')
end

if max(p.MaskChannelIndex) > nChan || min(p.MaskChannelIndex)<1 ...
        || ~isequal(round(p.MaskChannelIndex),p.MaskChannelIndex)
    error('Invalid mask channel numbers specified! Check MaskChannelIndex input!!')
end

if length(p.MaskChannelIndex) ~= length(p.ChannelIndex)    
    error('MaskChannelIndedx and ChannelIndex are different lengths! You must specify a mask channel for every channel to background subtract!')    
end



%% ---------- Init ---------- %%

%Find the background mask creation process, and the shade correction
%process, as this should have already been performed.
iBMProc = find(cellfun(@(x)(isa(x,'BackgroundMasksProcess')),movieData.processes_),1);                          
iSCProc = find(cellfun(@(x)(isa(x,'ShadeCorrectionProcess')),movieData.processes_),1);                          

nChanCorr = length(p.ChannelIndex);

if isempty(iBMProc) || isempty(iSCProc)
    error('Background masking and shade correction have not yet been performed on this movie! Please run first!!')
else        
    %Check which channels have background masks, shade correction
    hasBM = cellfun(@(x)(~isempty(x)),movieData.processes_{iBMProc}.outFilePaths_(1,:));
    hasSC = cellfun(@(x)(~isempty(x)),movieData.processes_{iSCProc}.outFilePaths_(1,:));
    
    %Check that these are the same
    if ~all(hasSC(p.ChannelIndex))
        error('Each channel to be background subtracted must have shade correction! Please apply shade correction to all needed channels before running background subtraction!')
    end
    if ~all(hasBM(p.MaskChannelIndex))
        error('Each mask channel specified by MaskChannelIndex must have background masks! Create background masks prior to background subtraction!')
    end   
    
    %Set the input directories to be the output directories of the
    %shade correction, and the correction images are the background masks
    %from the channels specified by MaskChannelIndex.
    for j = 1:nChanCorr
        
        movieData.processes_{iProc}.setCorrectionImagePath(p.ChannelIndex(j),...
                movieData.processes_{iBMProc}.outFilePaths_{p.MaskChannelIndex(j)});               
        
        movieData.processes_{iProc}.setInImagePath(p.ChannelIndex(j),...
            movieData.processes_{iSCProc}.outFilePaths_{1,p.ChannelIndex(j)});    
        
    end
end

%Set up the directories for corrected images as sub-directories of the
%output directory, and specify the directories for the images to be
%corrected in the movieData.


outFilePaths = cell(2,nChan);
for iChan = p.ChannelIndex;
    
    %Create string for current directory
    outFilePaths{1,iChan} = [p.OutputDirectory filesep dName num2str(iChan)];    
    
    %Check/create directory (checking avoids warning about existing
    %directory)
    mkClrDir(outFilePaths{1,iChan});
       
    outFilePaths{2,iChan} = [p.OutputDirectory filesep saveName num2str(iChan) '.mat'];   
end

%Save this in the process object
movieData.processes_{iProc}.setOutFilePaths(outFilePaths);

%Get image file names for input images (shade corrected images)
inNames =  movieData.processes_{iProc}.getInImageFileNames(p.ChannelIndex);

nImages = movieData.nFrames_;
nImTot = nImages * nChanCorr;



%% ---- Background Subtraction ---- %%
%Go through each image and subtract the average value behind the background
%mask from the image.


disp('Starting background subtraction...')

if ~p.BatchMode
    wtBar = waitbar(0,['Please wait, background subtracting channel ' num2str(p.ChannelIndex(1)) ' ...']);        
end        

backgroundValues = cell(1,nChan);


%Go through each image and apply the appropriate shade correction
for iChan = 1:nChanCorr
    
    inDir  = movieData.processes_{iProc}.inFilePaths_{1,p.ChannelIndex(iChan)};    
    outDir = movieData.processes_{iProc}.outFilePaths_{1,p.ChannelIndex(iChan)};    
    corrDir = movieData.processes_{iProc}.inFilePaths_{2,p.ChannelIndex(iChan)};
    
    bakNames = movieData.processes_{iBMProc}.getOutMaskFileNames(p.MaskChannelIndex(iChan));


    if ~p.BatchMode        
        waitbar((iChan-1)*nImages / nImTot,wtBar,['Please wait, background subtracting channel ' num2str(p.ChannelIndex(iChan)) ' ...']);        
    end        

    disp(['Background subtracting channel ' num2str(p.ChannelIndex(iChan)) '...'])
    disp(['Background subtracting images from "' inDir '", results will be stored in ' outDir]);     
    disp(['Using background masks from ' corrDir])
    
    
    backgroundValues{iChan} = zeros(1,nImages);
    
    for iImage = 1:nImages
    
        currIm = imread([inDir filesep inNames{iChan}{iImage}]);
        ogClass = class(currIm); %Determine class so it is not altered later
        
        %Load the background mask
        currBackMask = imread([corrDir filesep bakNames{1}{iImage}]);
        
        %Get average background intensity 
        backgroundValues{iChan}(iImage) = ...
            mean(double(currIm(currBackMask(:))));        
        currIm = double(currIm) - backgroundValues{iChan}(iImage);
        currIm(currIm < 0) = 0; %Clip negative values to zero
        currIm = cast(currIm,ogClass); %Return to original data type

        %Write it to disk        
        imwrite(currIm,[outDir filesep pString ...
            inNames{iChan}{iImage}]);                

        if ~p.BatchMode && mod(iImage,5)
            %Update the waitbar occasionally to minimize slowdown
            waitbar((iImage + (iChan-1)*nImages) / nImTot,wtBar)
        end                        
         

    end
end

if ~p.BatchMode && ishandle(wtBar)
    close(wtBar)
end


%% ------------- Output ------- %%

disp('Saving results...')

%Save the values that were subtracted from each frame.
for i = 1:nChanCorr
   subtractedValues = backgroundValues{i}; %#ok<NASGU>
   save(movieData.processes_{iProc}.outFilePaths_{2,p.ChannelIndex(i)},'subtractedValues');
end


%Log the correction in the movieData object and save it

movieData.processes_{iProc}.setDateTime;
movieData.save;

disp('Finished Correcting!')

