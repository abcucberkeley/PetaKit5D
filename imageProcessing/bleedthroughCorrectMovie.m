function movieData = bleedthroughCorrectMovie(movieData,varargin)
%BLEEDTHROUGHCORRECTMOVIE corrects for bleedthrough of other fluorophores in the input movie.
% 
% movieData = bleedthroughCorrectMovie(movieData)
% 
% movieData = bleedthroughCorrectMovie(movieData,paramsIn)
% 
%
% This function corrects for bleedthrough of other fluorophores into a
% particular channel in the input movie. This is done using bleedthrough
% coefficients calculated from a "bleedthrough movie" where only one
% fluorophore is present, using processBleedthroughMovie.m. These
% coefficients are then used, in combination with images of the fluorophore
% channels which are bleeding into the image to be corrected, to remove the
% effects of bleedthrough.
% 
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
%       sub-directory called "bleedthrough_corrected_images"
%
%       ('ChannelIndex'-> Positive integer scalar)
%       The integer index of the channel to perform bleedthrough/crosstalk
%       correction on. This index corresponds to the channel directories
%       location in the cell array movieData.channels). If not
%       input, the user will be asked to select from the movie's channels.       
%
%       ('Coefficients' -> n x 2 array where n is the number of channels)
%       The coefficients for each channel in the same order. The first
%       column gives the bleedthrough coefficients. The second channel
%       gives the cross-talk coefficients
%       This should be the average coefficient produced by 
%       calculateMovieBleedthrough.m
%
%       ('BatchMode' -> True/False)
%       If this option value is set to true, all graphical output is
%       suppressed.
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
% Hunter Elliott
% 11/2009
% Revamped 6/2010
%

%% ------ Parameters ------- %%

pString = 'btc_'; %The string to prepend before the bleedthrough-corrected image directory & channel name
dName = 'bleedthrough_corrected_images_for_channel_';%String for naming the directories for each corrected channel

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous bleedthrough correction process
iProc = movieData.getProcessIndex('BleedthroughCorrectionProcess',1);                          

%If the process doesn't exist, create it with default settings.
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(BleedthroughCorrectionProcess(movieData,movieData.outputDirectory_));                                                                                                 
end

%Parse input, store in parameter structure
p = parseProcessParams(movieData.processes_{iProc},paramsIn);

%% Initialization

%Ask the user for the image channel to correct if not input
if isempty(p.ChannelIndex)  
    if p.BatchMode
        error('In batch mode, you must specify the channel to perform bleedthrough correction on!')
    else
        %As the user to select a channel.
        p.ChannelIndex = selectMovieChannels(movieData,0,'Select the channel to bleedthrough correct:');        
    end    
elseif length(p.ChannelIndex) > 1
    % SB: this could be improved by storing a nx2xn matrix of coefficients
    % and  inverting this matrix for correcting multiple channels
    error('Only one channel may be bleedthrough/crosstalk corrected at a time!')
end

assert(movieData.processes_{iProc}.checkChanNum(p.ChannelIndex),...
    'Invalid channel number specified! Check ChannelIndex input!!')

%Ask the user for the channels which bleed into the images to be corrected, if not input
assert(~isempty(find(p.Coefficients,1)),'No bleedthrough coefficients');
bleedChannelIndex = find(p.Coefficients(:,1)~=0);
crossChannelIndex = find(p.Coefficients(:,2)~=0);
corrChannelIndex  = unique(vertcat(bleedChannelIndex,crossChannelIndex));

% If using the output of an existing process
if ~isempty(p.ProcessIndex)
    % Check which channels have been processed
    nChan = numel(movieData.channels_);
    nProc = numel(p.ProcessIndex);
    inChan = unique(vertcat(p.ChannelIndex,corrChannelIndex));
    
    hasOutput = false(nProc,nChan);    
    for i=1:nProc
        inProcIndex = p.ProcessIndex(i);        
        assert(isa(movieData.processes_{inProcIndex},'ImageProcessingProcess'));
                
        hasOutput(i,:)= movieData.processes_{inProcIndex}.checkChannelOutput;
    end
    assert(all(sum(hasOutput(:,inChan),1)),...
        'The channel to be corrected and the bleedthrough/crosstalk channels must all have been processed prior to bleedthrough correction!');    
    
    inProc = cell(1,nChan);
    for i=inChan', inProc{i} = movieData.processes_{p.ProcessIndex(find(hasOutput(:,i),1,'last'))}; end
end


%% Bleedthrough correction
disp('Starting bleedthrough/crosstalk correction...')

movieData.processes_{iProc}.setInFilePaths(cell(2,nChan));
%Retrieve the paths and names of the input images   
if isempty(p.ProcessIndex)
    imPaths = movieData.getChannelPaths;
    inNames = movieData.getImageFileNames(p.ChannelIndex);
else
    imPaths = cell(1,nChan);
    for i=inChan'
        imPaths{i} = inProc{i}.outFilePaths_{1, i};
    end
    inNames = inProc{p.ChannelIndex}.getOutImageFileNames(p.ChannelIndex);
end

%Log the input and bleed images in the process
movieData.processes_{iProc}.setInImagePath(p.ChannelIndex,imPaths{1,p.ChannelIndex});
movieData.processes_{iProc}.setCorrectionImagePath(corrChannelIndex,imPaths(1,corrChannelIndex));

% Set the output directory
outDir = [p.OutputDirectory filesep dName num2str(p.ChannelIndex)]; %Corrected images
mkClrDir(outDir);
movieData.processes_{iProc}.setOutImagePath(p.ChannelIndex,outDir);


%% -------------- Apply bleedthrough correction ------------%%
%Applies the bleedthrough correction from above to each selected channel

disp('Applying bleedthrough/crosstalk correction to images...')

%Go through each image and apply the appropriate bleedthrough correction
if ~p.BatchMode
    wtBar = waitbar(0,['Please wait, bleedthrough/crosstalk correcting channel ' num2str(p.ChannelIndex(1)) ' ...']);        
end        

nImages = movieData.nFrames_;

disp(['Correcting channel ' num2str(p.ChannelIndex) '...']);
disp(['Correcting images from ' imPaths{p.ChannelIndex}]);
disp(['Storing results in ' outDir]);

corrCoefficients = vertcat(p.Coefficients(bleedChannelIndex,1),p.Coefficients(crossChannelIndex,2));

for j = 1:length(corrChannelIndex)
    if j<=numel(bleedChannelIndex), type = 'bleedthrough'; else type='crosstalk'; end
    disp(['Correcting ' type ' from channel ' num2str(corrChannelIndex(j)) ','])
    disp(['using images from ' imPaths{corrChannelIndex(j)}])        
    disp(['...using ' type ' coefficient of ' num2str(corrCoefficients(j))])
end

for iImage = 1:nImages
    
    
    %Load the image to be corrected
    if isempty(p.ProcessIndex)
        currIm=movieData.channels_(p.ChannelIndex).loadImage(iImage);
    else
        currIm=inProc{p.ChannelIndex}.loadChannelOutput(p.ChannelIndex,iImage);
    end
    
    %Check the bit-depth of the image
    ogClass = class(currIm);
    currIm = double(currIm);

    for iCorr = 1:numel(corrChannelIndex)        
        
        %Load the correction image
        if isempty(p.ProcessIndex)
            currCorrIm=double(movieData.channels_(corrChannelIndex(iCorr)).loadImage(iImage));
        else
            currCorrIm=double(inProc{corrChannelIndex(iCorr)}.loadChannelOutput(corrChannelIndex(iCorr),iImage));
        end
%         currBleedIm = double(imread([bleedImDir{iBleed} filesep bleedImNames{iBleed}{iImage}]));

        %Subtract the bleedthrough from this channel
        currIm = currIm - (currCorrIm * corrCoefficients(iCorr));
        
        %Remove negative values (these usually occur in the background)
        currIm(currIm < 0) = 0;
        
        if ~any(currIm > 0)
            error('Please check the specified bleedthrough coefficients: The specified correction results in completely blank images!')
        end
        
    end  


    %Cast to original class
    currIm = cast(currIm,ogClass);

    %Write it to disk    
    imwrite(currIm,[outDir filesep pString inNames{1}{iImage}]);

    if ~p.BatchMode && mod(iImage,5)
        %Update the waitbar occasionally to minimize slowdown
        waitbar(iImage / nImages,wtBar)
    end                        


end

if ~p.BatchMode && ishandle(wtBar)
    close(wtBar)
end


%% ------------- Output ------- %%

disp('Saving results...')

%Log the correction in the movieData object and save it
movieData.processes_{iProc}.setDateTime;
movieData.save;

disp('Finished correcting bleedthrough/crosstalk!')