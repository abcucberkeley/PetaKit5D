% 
% The following script tests the steerableDetector3D on a synthetic 3D volume 
% with lines drawn in random orientations
% 
% Author: Deepak Roy Chittajallu (Created on Aug 31, 2012)
%

clc
clear
close all

%*************************************************************************
%                               PARAMETERS
%*************************************************************************

    numLines = 100;
    
    meanLineWidth = 2;
    stdLineWidth = 0.5;
    
%*************************************************************************

im = zeros( 200, 200, 30 );

% normal distributons for foreground and background pixels
fgMeanVar = [ 200, 20 ];
bgMeanVar = [ 180, 20 ];
fgGmObj = gmdistribution( fgMeanVar(:,1), reshape(fgMeanVar(:,2), [1,1,size(fgMeanVar,1)]) );
bgGmObj = gmdistribution( bgMeanVar(:,1), reshape(bgMeanVar(:,2), [1,1,size(bgMeanVar,1)]) );
im(:) = random(bgGmObj,numel(im));

% generate lines randomly
fprintf( '\nGenerating %d Lines in Random Orientations: \n', numLines );

imsize = size( im );
[X,Y,Z] = meshgrid(1:imsize(2), 1:imsize(1), 1:imsize(3));

try
    ptIm = gpuArray( [ X(:), Y(:), Z(:) ] );
catch err
    ptIm = [ X(:), Y(:), Z(:) ]; % you probably dont have a gpu, the code will be slow
end

lineMask = false( imsize );
lineWidthVec = zeros( numLines, 1 );

for i = 1:numLines   
    
    fprintf( ' %.3d ', i );
    
    if mod( i, 15 ) == 0
        fprintf( '\n' );
    end
    
    % generate a random point in volume
    ptRandInd = floor( rand * numel(im) );
    ptRef = ptIm( ptRandInd, : );
    
    % generate random orientation vector
    thetax = rand * pi;
    thetaxy = (0.05 + 0.05 * randn ) * pi;
    v = [ cos(thetaxy) * cos(thetax), cos(thetaxy) * sin(thetax), sin(thetaxy) ];
    
    % generate random line width
    randLineWidth = meanLineWidth + stdLineWidth * randn;
    
    if randLineWidth <= 0 
        randLineWidth = meanLineWidth;
    end
    
    lineWidthVec(i) = randLineWidth;
    
    % generate line
    refvec = ptIm - repmat( ptRef, numel(im), 1 );
    perp = refvec - (refvec * v') * v;
    sqDist = sum( perp .* perp, 2 );
    lineMask( sqDist <= randLineWidth^2 ) = true;
    
end

fprintf( '\n' );

im( lineMask ) = random(fgGmObj, numel( find( lineMask ) ));

% Run steerable detector to enhance the lines
fprintf( '\nRunning steerable detector at multiple scales on %d x %d x %d sized volume ...\n', imsize(2), imsize(1), imsize(3) );
sigmaTrialValues = meanLineWidth + stdLineWidth * (-2:1:2);
sigmaTrialValues( sigmaTrialValues <= 0 ) = [];

for i = 1:numel( sigmaTrialValues )
    
    fprintf( '\n\t%d/%d: Trying sigma value of %.2f ... ', i, numel( sigmaTrialValues ), sigmaTrialValues(i) );   
    
    tic
    [curRes, cutTheta, curNms] = steerableDetector3D(im, 1, sigmaTrialValues(i));
    timeElapsed = toc;
    
    fprintf( 'It took %.2f seconds\n', timeElapsed );   
    
    %curRes = sigmaTrialValues(i) * curRes; % scale normalization
    
    if i == 1        
        res = curRes;
        nms = curNms;
        theta = cutTheta;
        pixelScaleMap = ones( size(res) );
    else
        indBetter = curRes > res;
        res(indBetter) = curRes(indBetter);
        nms(indBetter) = curNms(indBetter);
        theta.x1(indBetter) = cutTheta.x1(indBetter);
        theta.x2(indBetter) = cutTheta.x2(indBetter);
        theta.x3(indBetter) = cutTheta.x3(indBetter);
        pixelScaleMap(indBetter) = i;
    end
end

imLineSegMask = res > thresholdOtsu(res);
imLineSegRGBMask = zeros( [size(imLineSegMask), 3] );
imLineSegRGBMask(:,:,:,1) = imLineSegMask;

pixelScaleMap( ~imLineSegMask ) = 0;
imScaleRGB = reshape( label2rgb( pixelScaleMap(:), 'jet', 'k' ) / 255.0, size(imLineSegRGBMask) );

% display
imseriesshow( im );
set( gcf, 'Name', 'Image with Randomly Generated Lines' );

imseriesmaskshowrgb( res, {imLineSegRGBMask, imScaleRGB} );
colorbar;
set( gcf, 'Name', 'Response of Steerable Detector with Line Mask and Scale Map' );

imseriesmaskshowrgb( nms, {imLineSegRGBMask, imScaleRGB} );
colorbar;
set( gcf, 'Name', 'Result of Non-maximal suppression with Line Mask and Scale Map' );
