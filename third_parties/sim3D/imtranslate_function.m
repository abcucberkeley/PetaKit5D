function [B,RB] = imtranslate_function(varargin)
%IMTRANSLATE Translate image.
%   B = IMTRANSLATE(A,TRANSLATION) translates image A by a translation
%   vector TRANSLATION. TRANSLATION is of the form [TX TY] for 2-D inputs,
%   and [TX TY TZ] for 3-D inputs. If TRANSLATION is a two-element vector
%   and A has more than two dimensions, a 2-D translation is applied to A
%   one plane at a time. TRANSLATION can be fractional.
%
%   [B,RB] = IMTRANSLATE(A,RA,TRANSLATION) translates the spatially
%   referenced image defined by A and RA by a translation vector
%   TRANSLATION. TRANSLATION is in the world coordinate system. The output
%   is a translated spatially referenced image defined by B and RB.
%
%   B = IMTRANSLATE(A,TRANSLATION,METHOD) translates image A, using the interpolation
%   method specified by METHOD. METHOD is a string that can have one of the
%   following values.
%
%        'nearest'      Nearest neighbor interpolation
%
%        'linear'       Linear interpolation
%
%        'cubic'        Cubic interpolation. Note: This interpolation
%                       method can produce pixel values outside the original
%                       range.
%    Default: The default interpolation is 'linear' for numeric inputs and
%    'nearest' for categorical inputs.
%    Categorical inputs only support the 'nearest' interpolation type.
%
%   [___] = IMTRANSLATE(___,Name, Value,...) translates the input image
%   using name-value pairs to control various aspects of the translation.
%
%   Parameters include:
%
%      'OutputView' - A string that defines the output world limits. Valid
%                     strings are 'same' and 'full'. When 'OutputView' is
%                     'same', the output world limits are the same as the
%                     input image. When 'OutputView' is 'full', the world
%                     limits are the bounding rectangle that includes both the
%                     input image and the translated output image.
%
%                     Default value: 'same'
%
%      'FillValues'   An array containing one or several fill values.
%                     Fill values are used for output pixels when the
%                     corresponding inverse transformed location in the
%                     input image is completely outside the input image
%                     boundaries.
%
%                     If A is 2-D then 'FillValues' must be a scalar. If A
%                     is 3-D and TRANSLATION is a 3-element vector, then
%                     'FillValues' must be a scalar. If A is N-D and
%                     TRANSLATION is a 2-element vector, then 'FillValues'
%                     may be either scalar or an array whose size matches
%                     dimensions 3 to N of A. For example, if A is a uint8
%                     RGB image that is 200-by-200-by-3, then 'FillValues'
%                     can be a scalar or a 3-by-1 array. In this RGB image
%                     example, possibilities for 'FillValues' include:
%
%                           0                 - fill with black
%                           [0;0;0]           - also fill with black
%                           255               - fill with white
%                           [255;255;255]     - also fill with white
%                           [0;0;255]         - fill with blue
%                           [255;255;0]       - fill with yellow
%
%                     If A is 4-D with size 200-by-200-by-3-by-10, then
%                     FillValues' can be a scalar or a 3-by-10 array.
%
%                     When A is categorical input, FillValues can only be a
%                     scalar and can take one of the following values:
%                        - Valid category in input data specified as string
%                        or character array.
%
%                        - missing, which corresponds  to <undefined>
%                        category (default)
%
%                     Default: 0 for numeric and missing(<undefined>)
%                     for categorical
%
%    Class Support
%    -------------
%    A can be of any nonsparse, numeric (except uint64 and int64) or
%    categorical. A can also be logical. TRANSLATION is a nonsparse, real
%    valued numeric vector.  The class of B is the same as the class of A.
%    RA and RB are spatial referencing objects of class imref2d or imref3d.
%
%    Notes
%    -----
%
%    1) IMTRANSLATE is optimized for integrally valued TRANSLATION vectors.
%
%    2) When 'OutputView' is 'full' and translation is a fractional number
%    of pixels, the world limits of the output spatial referencing object
%    Rout are expanded to the nearest full pixel increment such that Rout
%    contains both the original and translated images at the same
%    resolution as the input image A. The additional image extent in each
%    is added on one side of the image, in the direction that the
%    translation vector points. For example, when translation is
%    fractional and positive in both X and Y, then the maximum of
%    XWorldLimits and YWorldLimits is expanded to enclose the 'full'
%    bounding rectangle at the resolution of the input image.
%
%   Example 1
%   ---------
%   % Translate image I by 5.3 pixels in X and -10.1 pixels in Y.
%
%       I = imread('pout.tif');
%       J = imtranslate(I,[5.3, -10.1],'FillValues',255);
%       figure, montage({I, J});
%
%   Example 2
%   ---------
%   % Translate a 3-D MRI dataset. Use 'OutputView' to obtain full
%   % translated image volume without clipping.
%
%       s = load('mri');
%       mriVolume = squeeze(s.D);
%       sizeIn = size(mriVolume);
%       hFigOriginal = figure;
%       hAxOriginal  = axes;
%       slice(double(mriVolume),sizeIn(2)/2,sizeIn(1)/2,sizeIn(3)/2);
%       grid on, shading interp, colormap gray
%
%       % Apply a translation in the X,Y direction
%       mriVolumeTranslated = imtranslate(mriVolume,[40,30,0],'OutputView','full');
%
%       % Visualize axial slice plane taken through center of volume
%       sliceIndex = round(sizeIn(3)/2);
%       axialSliceOriginal   = mriVolume(:,:,sliceIndex);
%       axialSliceTranslated = mriVolumeTranslated(:,:,sliceIndex);
%
%       imshowpair(axialSliceOriginal,axialSliceTranslated,'montage');
%
%   See also IMRESIZE, IMROTATE, IMWARP.

%   Copyright 2013-2019 The MathWorks, Inc.

[R_A, varargin] = preparseSpatialReferencingObjects(varargin{:});

[A,translation,method,outputView,fillValues,catConverter, isInputCategorical] = parseInputs(varargin{:});

inputSpatialReferencingNotSpecified = isempty(R_A);
if inputSpatialReferencingNotSpecified
    if isa(R_A,'imref3d')
        R_A = imref3d(size(A));
    else
        R_A = imref2d(size(A));
    end
else
    % Check agreement of input spatial referencing object with input image.
    checkSpatialRefAgreementWithInputImage(A,R_A);
end

is2DProblem = isequal(numel(translation),2);

integrallyValuedTranslation = inputSpatialReferencingNotSpecified && all(mod(translation,1) == 0);

if integrallyValuedTranslation && isreal(A)
    % As a performance optimization, we treat non-spatially referenced,
    % real valued problems with integral translations in all dimensions as
    % a special case.
    if is2DProblem
        [B,RB] = translateIntegerShift2D(A,R_A,translation,method,outputView,fillValues);
    else
        [B,RB] = translateIntegerShift3D(A,R_A,translation,method,outputView,fillValues);
    end
    
else
    
    if is2DProblem
        [B,RB] = translate2D(A,R_A,translation,method,outputView,fillValues);
    else
        [B,RB] = translate3D(A,R_A,translation,method,outputView,fillValues);
    end
    
end

if isInputCategorical
    B = catConverter.numeric2Categorical(B);
end

end

function [out,Rout] = translateIntegerShift2D(A,RA,translation,~,outputView,fillValues)

% This code path is a special case for non-spatially referenced cases in
% which the translation is integrally valued in all dimensions. We can
% avoid the computational cost of interpolation in these cases and form the
% output image with simple indexing.

Rout = computeOutputSpatialRef(RA,translation,outputView);

% Determine size of output image
inputSize = size(A);
if length(inputSize) < 3
    outputSize = Rout.ImageSize;
    numPlanes = 1;
else
    outputSize = [Rout.ImageSize, inputSize(3:end)];
    numPlanes = prod(inputSize(3:end));
end

% Pre-allocate output image to FillValue.
fillValueCastToOutputType = cast(fillValues,class(A));
if isscalar(fillValues)
    % This pre-allocation has to work with logical values as well as
    % numeric types.
    out = zeros(outputSize,'like',A);
    out(:) = fillValueCastToOutputType;
else
    out = zeros(outputSize,'like',A);
    for i = 1:length(fillValues)
        out(:,:,i) = fillValueCastToOutputType(i);
    end
end

[XWorldBoundingSubscripts,YWorldBoundingSubscripts] = Rout.intrinsicToWorld([1 Rout.ImageSize(2)],...
    [1 Rout.ImageSize(1)]);

UWorld = XWorldBoundingSubscripts - translation(1);
VWorld = YWorldBoundingSubscripts - translation(2);

% Figure out whether bounding rectangle of reverse mapped pixel centers
% includes any in bounds locations in the source image.
intrinsicSourceBoundingRectangleU = [UWorld(1) UWorld(2) UWorld(2) UWorld(1)];
intrinsicSourceBoundingRectangleV = [VWorld(1) VWorld(1) VWorld(2) VWorld(2)];
% Contains is true inside the world limits and the boundary of the world
% limits.
locationsInSourceMapToDestination = any(RA.contains(intrinsicSourceBoundingRectangleU,...
    intrinsicSourceBoundingRectangleV));

% If there are locations in the source image that map into the destination,
% use indexing to form the output image. Otherwise, return all fill values.
if locationsInSourceMapToDestination
    
    % Clip reverse mapped boundaries to boundaries that live entirely
    % within A.
    UWorldClippedToBounds = [max(1,UWorld(1)), min(RA.ImageSize(2),UWorld(2))];
    VWorldClippedToBounds = [max(1,VWorld(1)), min(RA.ImageSize(1),VWorld(2))];
    
    % At this point we know the locations in source that map into valid
    % locations in the destination image. We want to forward map these into
    % corresponding subscripts in our output image.
    NonFillOutputLocX = UWorldClippedToBounds + translation(1);computeOutputSpatialRef
    NonFillOutputLocY = VWorldClippedToBounds + translation(2);
    [outputR,outputC] = Rout.worldToSubscript(NonFillOutputLocX,NonFillOutputLocY);
    
    % Where the output locations map into valid, in-bounds locations in A,
    % assign the output values by simple indexing. No interpolation is
    % required since the translation is integrally valued.
    for i = 1:numPlanes
        out(outputR(1):outputR(2),outputC(1):outputC(2),i) = A(VWorldClippedToBounds(1):VWorldClippedToBounds(2),...
            UWorldClippedToBounds(1):UWorldClippedToBounds(2),i);
    end
    
end


end

function [out,Rout] = translateIntegerShift3D(A,RA,translation,~,outputView,fillValues)

% This code path is a special case for non-spatially referenced cases in
% which the translation is integrally valued in all dimensions. We can
% avoid the computational cost of interpolation in these cases and form the
% output image with simple indexing.

Rout = computeOutputSpatialRef(RA,translation,outputView);

if isgpuarray(A)
    fillValueCastToOutputType = gpuArray(fillValues);
else
    fillValueCastToOutputType = cast(fillValues,class(A));
end

% This pre-allocation has to work with logical values as well as
% numeric types.
out = zeros(Rout.ImageSize,'like',A);
out(:) = fillValueCastToOutputType;

[XWorldBoundingSubscripts,YWorldBoundingSubscripts,ZWorldBoundingSubscripts] = Rout.intrinsicToWorld([1 Rout.ImageSize(2)],...
    [1 Rout.ImageSize(1)],[1 Rout.ImageSize(3)]);

UWorld = XWorldBoundingSubscripts - translation(1);
VWorld = YWorldBoundingSubscripts - translation(2);
WWorld = ZWorldBoundingSubscripts - translation(3);

% Figure out whether bounding rectangle of reverse mapped pixel centers
% includes any in bounds locations in the source image. In specific the
% bounding cube, wind CCW in the UV plane on the lower W face of the cube,
% then repeat winding CCW in the UV plane on the upper W face of the cube.
intrinsicSourceBoundingCubeU = [UWorld(1), UWorld(2), UWorld(2), UWorld(1),...
                                UWorld(1), UWorld(2), UWorld(2), UWorld(1)];
                            
intrinsicSourceBoundingCubeV = [VWorld(1), VWorld(1), VWorld(2), VWorld(2),...
                                VWorld(1), VWorld(1), VWorld(2), VWorld(2)];
                            
intrinsicSourceBoundingCubeW = [WWorld(1), WWorld(1), WWorld(1), WWorld(1),...
                                WWorld(2), WWorld(2), WWorld(2), WWorld(2)];
                            
% Contains is true inside the world limits and the boundary of the world
% limits.
locationsInSourceMapToDestination = any(RA.contains(intrinsicSourceBoundingCubeU,...
                                                    intrinsicSourceBoundingCubeV,...
                                                    intrinsicSourceBoundingCubeW));
                                                    
% If there are locations in the source image that map into the destination,
% use indexing to form the output image. Otherwise, return all fill values.
if locationsInSourceMapToDestination
    
    % Clip reverse mapped boundaries to boundaries that live entirely
    % within A.
    UWorldClippedToBounds = [max(1,UWorld(1)), min(RA.ImageSize(2),UWorld(2))];
    VWorldClippedToBounds = [max(1,VWorld(1)), min(RA.ImageSize(1),VWorld(2))];
    WWorldClippedToBounds = [max(1,WWorld(1)), min(RA.ImageSize(3),WWorld(2))];

    % At this point we know the locations in source that map into valid
    % locations in the destination image. We want to forward map these into
    % corresponding subscripts in our output image.
    NonFillOutputLocX = UWorldClippedToBounds + translation(1);
    NonFillOutputLocY = VWorldClippedToBounds + translation(2);
    NonFillOutputLocZ = WWorldClippedToBounds + translation(3);

    [outputR,outputC,outputP] = Rout.worldToSubscript(NonFillOutputLocX,NonFillOutputLocY,NonFillOutputLocZ);
    
    % Where the output locations map into valid, in-bounds locations in A,
    % assign the output values by simple indexing. No interpolation is
    % required since the translation is integrally valued.
    out(outputR(1):outputR(2),outputC(1):outputC(2),outputP(1):outputP(2)) = A(VWorldClippedToBounds(1):VWorldClippedToBounds(2),...
        UWorldClippedToBounds(1):UWorldClippedToBounds(2),WWorldClippedToBounds(1):WWorldClippedToBounds(2));
    
end

end

function [out,Rout] = translate2D(A,RA,translation,method,outputView,fillValues)

Rout = computeOutputSpatialRef(RA,translation,outputView);

% Compute spatially referenced case as a general 2-D affine transformation.
tform = affine2d([1 0 0; 0 1 0; translation(1:2) 1]);
[out,Rout] = imwarp(A,RA,tform,method,'OutputView',Rout,'fillValues',fillValues,...
    'SmoothEdges', true);

end


function [out,Rout] = translate3D(A,RA,translation,method,outputView,fillValues)

Rout = computeOutputSpatialRef(RA,translation,outputView);

% Compute spatially referenced case as a general 3-D affine transformation.
tform = affine3d([1 0 0 0; 0 1 0 0; 0 0 1 0; translation(1:3) 1]);
[out,Rout] = imwarp(A,RA,tform,method,'OutputView',Rout,'fillValues',fillValues,...
    'SmoothEdges', false);

end


function Rout = computeOutputSpatialRef(RA,translation,outputView)

if strcmp(outputView,'same')
    Rout = RA;
else
    
    % imtranslate(___,'OutputView','full');
    [XWorldLimitsOut,numColsOutput] = computeFullExtentAndGridSizePerDimension(RA.XWorldLimits,...
        RA.PixelExtentInWorldX,translation(1));
    
    [YWorldLimitsOut,numRowsOutput] = computeFullExtentAndGridSizePerDimension(RA.YWorldLimits,...
        RA.PixelExtentInWorldY,translation(2));
    
    if ~isa(RA,'imref3d')
        
        Rout = imref2d([numRowsOutput numColsOutput],XWorldLimitsOut,YWorldLimitsOut);
             
    else
           
        [ZWorldLimitsOut,numPlanesOutput] = computeFullExtentAndGridSizePerDimension(RA.ZWorldLimits,...
            RA.PixelExtentInWorldZ,translation(3));
        
        Rout = imref3d([numRowsOutput numColsOutput numPlanesOutput],XWorldLimitsOut,YWorldLimitsOut,ZWorldLimitsOut);        
    end
    
end

end

function [worldLimitsOut,numPixelsInDimOutput] = computeFullExtentAndGridSizePerDimension(inputWorldLimits,inputWorldPixelExtentInDim,translationInDim)

% The full bounding rectangle is the bounding rectangle that
% includes the original and translated images
worldLimitsTranslated = inputWorldLimits+translationInDim;
minInDim = min(inputWorldLimits(1),worldLimitsTranslated(1));
maxInDim = max(inputWorldLimits(2),worldLimitsTranslated(2));

worldLimitsFullIdeal = [minInDim maxInDim];
idealFullExtent = diff(worldLimitsFullIdeal);

% Compute the number of pixels necessary to capture the entire full
% bounding box at the input image resolution. If the full extent is
% not evenly divisible by the input image resolution, use ceil to
% guarantee that we completely capture the full bounding box at the
% input image resolution.
numPixelsInDimOutput = ceil(idealFullExtent ./ inputWorldPixelExtentInDim);

% Compute the extent in world units of the output image, determined
% by the input image resolution and the number of pixels in the output
% image along each dimension.
outputImageExtentInDim = numPixelsInDimOutput*inputWorldPixelExtentInDim;

% If the ideal full image extent is not evenly divisible by the
% input image resolution, then the ceil will have added additional
% image extent. Compute the additional image extent.
addedImageExtentInDim  = outputImageExtentInDim - idealFullExtent;

% Add the additional image extent in each dimension on one side of
% the output image. Increase the full bounding box on the side that
% the translation vector points toward. This allows for a gradual
% transition as the translated image moves in sub-pixel increments.
if translationInDim >=0
    worldLimitsOut = worldLimitsFullIdeal + [0 addedImageExtentInDim];
else
    worldLimitsOut = worldLimitsFullIdeal + [-addedImageExtentInDim 0];
end

end


function [A,translation,method,outputView,fillValues,catConverter, isInputCategorical] = parseInputs(varargin)

catConverter = [];
isInputCategorical = false;
varargin = matlab.images.internal.stringToChar(varargin);

p = inputParser();
p.addRequired('A');
p.addRequired('TRANSLATION');
p.addParameter('OutputView', 'same')

if iscategorical(varargin{1})
    catConverter = images.internal.utils.CategoricalConverter(categories(varargin{1}));
    p.addParameter('FillValues', missing)
    p.addOptional('METHOD', 'nearest');
else
    p.addParameter('FillValues', 0)
    p.addOptional('METHOD', 'linear');
end

p.parse(varargin{:});

A               = p.Results.A;
translation     = double(p.Results.TRANSLATION);
method          = p.Results.METHOD;
outputView      = p.Results.OutputView;
fillValues      = p.Results.FillValues;

if any(strcmp(method,{'linear','cubic'}))
    method = strcat('bi',method);
end

if iscategorical(A)
    % For categorical inputs,
    % 1. 'nearest' is the only interpolation method.
    % 2.  Default Fill values will be set to missing, which corresponds to
    % '<undefined>' label in the output categorical result. Any other
    % FillValue results in an error.
   
    % Get corresponding numeric value for the classname
    fillValues = catConverter.getNumericValue(fillValues);
    
    A = catConverter.categorical2Numeric(A);
    isInputCategorical = true;
end

end

function [R_A,varargin] = preparseSpatialReferencingObjects(varargin)
% This should be abstracted into a package

if (nargin > 1) && (isa(varargin{2},'imref2d') || isa(varargin{2},'imref3d'))
    validateattributes(varargin{2},{'imref2d','imref3d'},{'scalar','nonempty'},'imwarp','RA');
    R_A = varargin{2};
    varargin(2) = [];
else
    % We don't want to actually assign the default spatial referencing
    % object until the rest of the input arguments have been validated.
    % Assign empty spatial referencing arguments as a flag that we need to
    % assign the identity spatial referencing object after input
    % parsing/validation has finished.
    translation = varargin{2};
    if (numel(translation) == 2)
        R_A = imref2d.empty();
    else
        R_A = imref3d.empty();
    end
    
end

end





