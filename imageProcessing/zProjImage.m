function imageProj = zProjImage(imageFile,projTypes)

%% Input

%ask user for images
if nargin < 1 || isempty(imageFile)
    [fName,dirName] = uigetfile('*.tif','Specify image file');
    imageFile = fullfile(dirName,fName);
else
    if iscell(imageFile)
        [fpath,fname,fno,fext]=getFilenameBody(imageFile{1});
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    elseif ischar(imageFile)
        [fpath,fname,fno,fext]=getFilenameBody(imageFile);
        dirName=[fpath,filesep];
        fName=[fname,fno,fext];
    end
end

%if input is valid ...
if(isa(fName,'char') && isa(dirName,'char'))
    
    %get number of frames in stack
    numFrames = length(imfinfo(imageFile));
    
else %else, exit
    
    disp('--zProjImage: Bad file selection');
    return
    
end

if nargin < 2 || isempty(projTypes)
    projTypes = 'max';
end
numProj = length(projTypes);

%% Projection

%read first image to get image size and reserve memory
image1 = imread(imageFile,1);
imagesAll = repmat(image1,[1 1 numFrames]);
imageProj = double(repmat(image1,[1 1 numProj]));

%read all images
for iFrame = 2 : numFrames
    imagesAll(:,:,iFrame) = imread(imageFile,iFrame);
end
imagesAll = double(imagesAll);

%do different types of projections
for iProj = 1 : numProj
   
    %get current projection type
    projCurr = projTypes{iProj};
    
    %do projection
    if strcmp(projCurr,'max')
        
        imageProjTmp = max(imagesAll,[],3);
        imwrite(uint16(imageProjTmp),fullfile(dirName,['MAX_' fName]),'tif');
        
    elseif strcmp(projCurr,'ave')
        
        imageProjTmp = nanmean(imagesAll,3); 
        imwrite(uint16(imageProjTmp),fullfile(dirName,['AVE_' fName]),'tif');
        
    elseif strcmp(projCurr,'med')
        
        imageProjTmp = nanmedian(imagesAll,3);
        imwrite(uint16(imageProjTmp),fullfile(dirName,['MED_' fName]),'tif');
        
    end
    
    imageProj(:,:,iProj) = imageProjTmp;
    
end

%% ~~~ the end ~~~
