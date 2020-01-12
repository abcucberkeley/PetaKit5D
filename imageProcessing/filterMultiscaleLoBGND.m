function [ imMultiscaleLoBGResponse, varargout ] = filterMultiscaleLoBGND( imInput, sigmaValues, rho, varargin )
% filterMultiscaleLoBGND: An ND implementation of a multiscale Laplacian of Bi-gaussian (LoBG) Filter
% 
%     [ imMultiscaleLoGBResponse ] = filterMultiscaleLoGBND(imInput, sigmaValues, rho, varargin );
%     [ imMultiscaleLoGBResponse, pixelScaleMap ] = filterMultiscaleLoGBND(imInput, sigmaValues, rho, varargin );
% 
%     Required Input Arguments:
% 
%                            im: Input ND Image
%         
%                   sigmaValues: An array of scales (standard deviation of the 
%                                gaussian kernel) on which the LoBG filter 
%                                will be applied.
% 
%                                Caution: If the pixel spacing is not 
%                                specified then unit isotropic spacing is 
%                                assumed as a result of which the units of 
%                                sigmas will be assumed to be pixels. So if 
%                                you want to provide sigmas in physical
%                                space units then you should also provide
%                                the pixel spacing as an additional argument.
% 
%                           rho: ratio of sigmaBackground and sigma
% 
%     Optional Input Arguments:
% 
%                       spacing: pixel spacing of the input image.
%                                can be a scalar value or an array of size
%                                equal to the number of image dimensions.
%                                Default: 1 (isotropic spacing)                      
%                 
%               borderCondition: specifies the way in which the image is 
%                                padded at the borders. This argument is
%                                supplied as input to the functuin 'padarrayXT'. 
%                                Default: 'symmetric'
%                                Options: 'symmetric', 'replicate', 'circular', 
%                                         'antisymmetric', or a constant value
% 
%                     debugMode: true/false
%                                A bunch of stuff is printed in debug mode
% 
%                        UseGPU: true/false
%                                A flag that specifies whether or not to
%                                use the GPU for convolution.  
%                                   True - Uses GPU if installed
%                                   False - otherwise (default-value)
% 
%     Output Arguments:
% 
%       imMultiscaleLoGResponse: Response of the multiscale LoG Filter
% 
%                 pixelScaleMap: A map indicating for each pixel the
%                                scale/sigma at which the LoG response 
%                                was optimal accross the scale space
% 
%     Example: Comparison of LoG and LoBG responses on two adjacent
%     gaussian ridges of varying width
%         
%         g = @(x,sigma) ( exp(-x.^2/ (2*sigma^2)) );
%         ridge = @(x,w) ( g(x,w/2.5) ); 
%         multiridge = @(x,w) ( 2 * ridge(x,w) + 2 * ridge(x-1.5*w,2*w) + 0.1 * rand(1,numel(x)) );   
%         
%         w = 10; % width of ridge
%         sampleSpacing = w/100;
%         x = -8*w:sampleSpacing:8*w; 
%         rho = 0.1;
%         
%         % set sigma values: sigma = w/sqrt(12) is the optimal for ridge of width w
%         sigmaValues = (w/2) * 2.^[-1:1]; 
%         
%         figure; hold all;
%         strLegend = {};
% 
%         % plot ridge
%         plot( x, multiridge(x,w), '-' );
%         
%         strLegend{end+1} = 'Multiridge profile';
% 
%         % plot response of multiscale LoBG filter
%         [ multiscaleLoG, pixelScaleMap ] = filterMultiscaleLoGND( multiridge(x,w), sigmaValues, ...
%                                                                   'spacing', sampleSpacing );
%         plot( x, -multiscaleLoG, '-' ); 
%         strLegend{end+1} = 'Multiscale LoG Response';
% 
%         % plot response of multiscale LoBG filter
%         [ multiscaleLoBG, pixelScaleMap ] = filterMultiscaleLoBGND( multiridge(x,w), sigmaValues, rho, ...
%                                                                   'spacing', sampleSpacing );
%         plot( x, -multiscaleLoBG, '-' );          
%         strLegend{end+1} = 'Multiscale LoBG Response';
% 
%         title( 'Application of LoBG to two adjacent 1D Gaussian Ridges of varying widths', 'FontWeight', 'bold' );
%         legend( strLegend );      
% 
%   References:
% 
%   Xiao, C., M. Staring, et al. (2012). 
%   "A multiscale bi-Gaussian filter for adjacent curvilinear structures 
%   detection with application to vasculature images." 
%   IEEE Transactions on Image Processing, PP(99): 1-1.
% 
%   See: filerLoBGND.m
% 
%   Author: Deepak Roy Chittajallu
% 

    p = inputParser;
    p.CaseSensitive( false );
    p.addRequired( 'imInput', @(x) ( isnumeric(x) ) );    
    p.addRequired( 'sigmaValues', @(x) ( isnumeric(x) && numel(x) == max(size(x)) ) ); 
    p.addRequired( 'rho', @(x) ( isscalar(x) ) );
    p.parse( imInput, sigmaValues, rho );
    
    p.addParamValue( 'spacing', ones( 1, ndims(imInput) ), @(x) (isnumeric(x) && ismember(numel(x), [1,ndims(imInput)])) );
    p.addParamValue( 'borderCondition', 'symmetric', @(x) ( isscalar(x) || (ischar(x) && ismember( lower(x), { 'symmetric', 'replicate', 'antisymmetric' } ) ) ) );
    p.addParamValue( 'debugMode', false, @(x) (isscalar(x) && islogical(x)) );    
    p.addParamValue( 'UseGPU', false, @(x) (isscalar(x) && islogical(x)) );
    p.parse( imInput, sigmaValues, rho, varargin{:} );
    
    spacing = p.Results.spacing;
    borderCondition = p.Results.borderCondition;
    flagDebugMode = p.Results.debugMode;
    flagUseGPU = p.Results.UseGPU;
    
    
    % Run the LoG filter accross the scale space and record the optimal
    % response and the scale at which the optimal response was found for
    % each pixel
    if flagDebugMode 
       fprintf( '\nRunning LoBG filter at multiple scales on an image of size [ %s ] ...\n', ... 
                sprintf( ' %d ', size(imInput) ) );  
    end
        
    for i = 1:numel( sigmaValues )

        if flagDebugMode
            fprintf( '\n\t%d/%d: Trying sigma value of %.2f ... ', i, numel( sigmaValues ), sigmaValues(i) );   
            tic
        end

        [ imCurLoBGResponse ] = filterLoBGND( imInput, sigmaValues(i), rho, ... 
                                              'spacing', spacing, ...
                                              'borderCondition', borderCondition, ...
                                              'UseNormalizedDerivatives', true, ...
                                              'UseGPU', flagUseGPU );

        if flagDebugMode
            timeElapsed = toc;
            fprintf( 'It took %.2f seconds\n', timeElapsed );           
        end

        if i == 1

           imMultiscaleLoBGResponse = imCurLoBGResponse;
           pixelScaleMap = ones( size( imInput ) );

        else

            imBetterMask = imCurLoBGResponse < imMultiscaleLoBGResponse;
            imMultiscaleLoBGResponse( imBetterMask ) = imCurLoBGResponse( imBetterMask );
            pixelScaleMap( imBetterMask ) = i;

        end

    end
    
    if nargout > 1 
        varargout{1} = pixelScaleMap;
    end    
    
end
    
