function correctedImage = photobleachCorrection(images,masks)
% This function performs photobleach correction on a series of fluorescence
% images
%
%USAGE:
%       correctedImage = photobleachCorrection(images,masks)
%
%Inputs:
%       images - 3D matrix with the images to be corrected (x,y,#ofFrames)
%       masks  - same size as images
%
%Output:
%
%       correctedImages
%
%Marco Vilela, 2014
if size(images) ~= size(masks)
    error('Images and masks have different sizes')
end

[xL,yL,nFrame] = size(images);

SS  = ones(yL,1,nFrame);
SST = permute(SS,[2 1 3]);

totalF = cell2mat(arrayfun(@(x) SST(:,:,x)*(images(:,:,x).*masks(:,:,x))*SS(:,:,x),1:nFrame,'Unif',0));
countF = cell2mat(arrayfun(@(x) SST(:,:,x)*masks(:,:,x)*SS(:,:,x),1:nFrame,'Unif',0));
meanF  = totalF./countF;
time   = (0:nFrame-1)';

[fittedCurve,GOF] = fit(time,meanF','exp2');

%Testing consistency: variance of the data explained by the model - 
%Reference: Handbook of time series analysis. Wienheim:Wiley; 2006. p. 438?60

if GOF.rsquare < 0.8
    warndlg('WARNING:Crappy fit. Model explained less than 80% of the data variance')
    plot(fittedCurve,time,meanF')
end

correctedImage = reshape(cell2mat(arrayfun(@(x) images(:,:,x)*(1/fittedCurve(x)),1:nFrame,'Unif',0)),xL,yL,nFrame);
