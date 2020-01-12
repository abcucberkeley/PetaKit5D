function [ imMultiscaleLoGResponse, varargout ] = filterMultiscaleLoGND( imInput, sigmaValues, varargin )
% filterMultiscaleLoGND: An ND implementation of a multiscale Laplacian of Gaussian (LoG) Filter
% 
%     [ imMultiscaleLoGResponse ] = filterMultiscaleLoGND(imInput, sigmaValues, varargin );
%     [ imMultiscaleLoGResponse, pixelScaleMap ] = filterMultiscaleLoGND(imInput, sigmaValues, varargin );
% 
%     Required Input Arguments:
% 
%                            im: Input ND Image
%         
%                   sigmaValues: An array of scales (standard deviation of the 
%                                gaussian kernel) on which the LoG filter 
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
%         UseNormalizedGaussian: true/false
%                                specifies whether the guassian kernel should
%                                be normalized or not
%                                Default: true
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
%     Example: Application of LoG to flat 1D Ridges of varying widths
%         
%         step = @(x) ( double( x >= 0 ) ); % 1-D step edge
%         ridge = @(x,w) ( step(x + w/2) - step(x - w/2) ); % 1-D ridge of width w
%         multiridge = @(x,w) ( ridge(x+2.75*w,w/2) + ridge(x,w) + ridge(x-4*w,2*w) );   
%         
%         w = 10; % width of ridge
%         sampleSpacing = w/100;
%         x = -8*w:sampleSpacing:8*w; 
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
%         % plot response of multiscale LoG filter
%         [ multiscaleLoG, pixelScaleMap ] = filterMultiscaleLoGND( multiridge(x,w), sigmaValues, ...
%                                                                   'spacing', sampleSpacing );
%         plot( x, -multiscaleLoG, '-' ); 
%         
%         strLegend{end+1} = 'Multiscale LoG Response';
% 
%         title( 'Application of LoG to Flat 1D Ridges of varying widths', 'FontWeight', 'bold' );
%         legend( strLegend );        
% 
% 
%     Author: Deepak Roy Chittajallu
% 

    p = inputParser;
    p.CaseSensitive( false );
    p.addRequired( 'imInput', @(x) ( isnumeric(x) ) );    
    p.addRequired( 'sigmaValues', @(x) ( isnumeric(x) && numel(x) == max(size(x)) ) ); 
    p.parse( imInput, sigmaValues );
    
    p.addParamValue( 'spacing', ones( 1, ndims(imInput) ), @(x) (isnumeric(x) && ismember(numel(x), [1,ndims(imInput)])) );
    p.addParamValue( 'borderCondition', 'symmetric', @(x) ( isscalar(x) || (ischar(x) && ismember( lower(x), { 'symmetric', 'replicate', 'antisymmetric' } ) ) ) );
    p.addParamValue( 'UseNormalizedGaussian', true, @(x) (isscalar(x) && islogical(x)) );    
    p.addParamValue( 'debugMode', false, @(x) (isscalar(x) && islogical(x)) );    
    p.addParamValue( 'UseGPU', false, @(x) (isscalar(x) && islogical(x)) );
    p.parse( imInput, sigmaValues, varargin{:} );
    
    spacing = p.Results.spacing;
    borderCondition = p.Results.borderCondition;
    flagNormalizeGaussian = p.Results.UseNormalizedGaussian;
    flagDebugMode = p.Results.debugMode;
    flagUseGPU = p.Results.UseGPU;
    
    % Run the LoG filter accross the scale space and record the optimal
    % response and the scale at which the optimal response was found for
    % each pixel
    if flagDebugMode 
       fprintf( '\nRunning LoG filter at multiple scales on an image of size [ %s ] ...\n', ... 
                sprintf( ' %d ', size(imInput) ) );  
    end
        
    for i = 1:numel( sigmaValues )

        if flagDebugMode
            fprintf( '\n\t%d/%d: Trying sigma value of %.2f for blobs of diameter %.2f ... ', ...
                     i, numel( sigmaValues ), sigmaValues(i), ...
                     sigmaValues(i) * 2 * sqrt(ndims(imInput)) );   
            tic
        end

        [ imCurLoGResponse ] = filterLoGND( imInput, sigmaValues(i), ... 
                                            'spacing', spacing, ...
                                            'borderCondition', borderCondition, ...
                                            'UseNormalizedDerivatives', true, ...
                                            'UseNormalizedGaussian', flagNormalizeGaussian, ...
                                            'UseGPU', flagUseGPU );

        if flagDebugMode
            timeElapsed = toc;
            fprintf( 'It took %.2f seconds\n', timeElapsed );           
        end

        if i == 1

           imMultiscaleLoGResponse = imCurLoGResponse;
           pixelScaleMap = ones( size( imInput ) );

        else

            imBetterMask = imCurLoGResponse < imMultiscaleLoGResponse;
            imMultiscaleLoGResponse( imBetterMask ) = imCurLoGResponse( imBetterMask );
            pixelScaleMap( imBetterMask ) = i;

        end

    end
    
    if nargout > 1 
        varargout{1} = pixelScaleMap;
    end
    
end

