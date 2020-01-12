function [ imLoBGResponse, varargout ] = filterLoBGND( imInput, sigma, rho, varargin )
% filterLoBGND: An ND implementation of the Laplacian of Bi-Gaussian Filter
% 
% This filter is shown to be more robust than the Laplacian of Gaussian 
% filter for detecting blobs under heavy overlap
% 
%     [ imLoBGResponse ] = filterLoBGND(im, sigma, varargin );
% 
%     Required Input Arguments:
% 
%                            im: Input ND Image
%         
%                         sigma: standard deviation of the gaussian
%                                should be a scalar value
% 
%                                Caution: If the pixel spacing is not 
%                                specified then unit isotropic spacing is 
%                                assumed as a result of which the units of 
%                                sigma will be assumed to be pixels. So if 
%                                you want to provide sigma in physical
%                                space units then you should also specify
%                                the pixel spacing.
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
%      UseNormalizedDerivatives: true/false
%                                specifies whether the gaussian derivatives
%                                should be scale-normalized.
%                                Default: false
%                                    
%         UseNormalizedGaussian: true/false
%                                specifies whether the guassian kernel should
%                                be normalized or not
%                                Default: true
%                                   
%                        UseGPU: true/false
%                                A flag that specifies whether or not to
%                                use the GPU for convolution. % 
%                                   True - Uses GPU if installed
%                                   False - otherwise (default-value)
%                                       
%     Output Arguments:
% 
%                         imLoBG: LoBG filtered image
% 
%                     lobgKernel: returns the Laplacian of Bigaussian kernel
%                                 (optional output argument)   
% 
% References:
% 
%   Xiao, C., M. Staring, et al. (2012). 
%   "A multiscale bi-Gaussian filter for adjacent curvilinear structures 
%    detection with application to vasculature images." 
%    IEEE Transactions on Image Processing, PP(99): 1-1.
% 
%     Example-1: Application of LoBG with rho=1 on a 1D ridge/blob (for rho=1 the 
%     response should be same as LoG filter)
%  
%         step = @(x) ( double( x >= 0 ) ); % 1-D step edge
%         ridge = @(x,w) ( step(x + w/2) - step(x - w/2) ); % 1-D ridge of width w
%         g = @(x,sigma) ( (sigma * sqrt(2*pi))^-1 * exp( -x.^2 / (2*sigma^2) ) ); % gaussian
%         dg = @(x,sigma) ( (-x/sigma^2) .* g(x,sigma) ); % first derivative of gaussian
%         d2g = @(x,sigma) ( ((x.^2 - sigma^2)/sigma^4) .* g(x,sigma) ); % second-derivative (laplacian) of gaussian
%         truelogresponse = @(x,w,sigma) ( dg(x + w/2,sigma) - dg(x - w/2, sigma) ); % true LoG response for ridge of width w
%         
%         w = 10; % width of ridge
%         sampleSpacing = w/100;
%         x = -4*w:sampleSpacing:4*w; 
%         datadims = numel( find( size(x) > 1 ) );
%         rho = 1.0;
% 
%         ridgesignal = ridge(x,w);
%         
%         % check the difference between the response of this filter and the true LoG response
%         figure; hold all;        
%         sigma = w/(2*sqrt(datadims)); % this gives maximal reponse at ridge-center
%         
%             % compute response of LoBG 
%             [lobgResponse, lobgKernel] = filterLoBGND( ridgesignal, sigma, rho, ...
%                                                        'spacing', sampleSpacing, ...
%                                                        'UseNormalizedDerivatives', true ); 
% 
%             % plot ridge
%             plot( x, ridgesignal, '-' );
% 
%             % plot scale-normalized (multiply by sigma^2) log kernel/filter
%             plot( x, -sigma^2 * d2g(x,sigma), '-' );
% 
%             % plot true scale-normalized (multiply by sigma^2) LoG response
%             plot( x, -sigma^2 * truelogresponse(x,w,sigma), '-' ); 
% 
%             % plot response of our LoBG filter
%             plot( x, -lobgResponse, '-' ); 
%         
%             title( sprintf( 'Comparison of our LoBG filter with rho=%.2f with the true LoG response for a 1D Ridge', rho ), 'FontWeight', 'bold' );
%             legend( sprintf( '1D Ridge (width = %d)', w ), ... 
%                     sprintf( 'True LoG Kernel (\\sigma = width/%.2f)', w/sigma), ...
%                     sprintf( 'LoBG Kernel (\\rho = %.2f)', rho), ...
%                     'True LoG Response', ...
%                     'LoBG Response' );          
%         
%     Example-2: Comparison of the Responses of LoG and LoBG on two adjacent 1D ridges
%         
%         step = @(x) ( double( x >= 0 ) ); % 1-D step edge
%         ridge = @(x,w) ( step(x + w/2) - step(x - w/2) ); % 1-D ridge of width w
% 
%         w = 10; % width of ridge
%         sampleSpacing = w/100;
%         x = -4*w:sampleSpacing:4*w; 
%         datadims = 1;
%         rho = 0.2;
% 
%         ridgesignal = ridge(x,w) + 2 * ridge(x-1.0*w,w);
%         sigma = w/(2*sqrt(datadims)); % this gives maximal reponse at ridge-center for LoG 
% 
%         % compute response of LoG 
%         [logResponse, logKernel] = filterLoGND( ridgesignal, sigma, ...
%                                                 'spacing', sampleSpacing, ...
%                                                 'UseNormalizedDerivatives', true );
% 
%         % compute response of LoBG 
%         [lobgResponse, lobgKernel] = filterLoBGND( ridgesignal, sigma, rho, ...
%                                                   'spacing', sampleSpacing, ...
%                                                   'UseNormalizedDerivatives', true ); 
%         % plot responses
%         figure; hold all;
% 
%         plot( x, ridgesignal, '-' ); 
%         plot( x, -logResponse, '-' );
%         plot( x, -lobgResponse, '-' );
% 
%         title( 'Comparison of LoG and LoBG reponses on two adjacent ridges', 'FontWeight', 'bold' );
%         legend(  'Adjacent Ridge Signal', ...
%                  sprintf('LoG Response (\\sigma = %.2f)', sigma), ...                     
%                  sprintf('LoBG Response (\\rho = %.2f)', rho) );
% 
%     Example-3: Application of LoBG to two adjacent 2D blobs
%  
%         blob = @(x,y,r) ( double(x.^2 + y.^2 - r.^2 <= 0) ); % 2D blob of radius r 
% 
%         r = 5; % radius of blob
%         sampleSpacing = r/100;
%         [x,y] = meshgrid(-4*r:sampleSpacing:4*r, -4*r:sampleSpacing:4*r); 
%         datadims = numel( find( size(x) > 1 ) );
%         rho = 0.2;
% 
%         blobsignal = double(blob(x-0.5*r,y,r) | blob(x+0.5*r,y,r));
%         
%         % check the difference between the response of this filter and the true LoG response
%         
%         % compute response of LoG 
%         [logResponse, logKernel] = filterLoGND( blobsignal, r/sqrt(datadims), ...
%                                                 'spacing', sampleSpacing * ones(1,2), ...
%                                                 'UseNormalizedDerivatives', true );
% 
%         % compute response of LoBG 
%         [lobgResponse, lobgKernel] = filterLoBGND( blobsignal, r, rho, ...
%                                                   'spacing', sampleSpacing * ones(1,2), ...
%                                                   'UseNormalizedDerivatives', true ); 
%         % plot responses
%         figure; 
% 
%         hprofile = floor(0.5*size(x,1));       
%         subplot(2,3,[1:3]); 
%         hold all;
%         plot( x(hprofile, :), blobsignal(hprofile, :), 'b-' ); 
%         plot( x(hprofile, :), -logResponse(hprofile, :), 'g-' );
%         plot( x(hprofile, :), -lobgResponse(hprofile, :), 'r-' );%           
%         hold off;
% 
%         title( 'Comparison of LoG and LoBG reponses on a 2D blob', 'FontWeight', 'bold' );
%         legend(  sprintf( 'Mid y-profile of 2D blob (r = %.2f)', r ), ...
%                  sprintf('LoG Response (\\sigma = %.2f)', sigma), ...                     
%                  sprintf('LoBG Response (\\rho = %.2f)', rho) );
%   
%         subplot(2,3,4);
%         imshow( blobsignal, [] );
%         hold on;
%         plot(1:size(x,2), hprofile*ones(1,size(x,2)), 'b-');
%         hold off;
%         title( sprintf( '2D blob with radius = %.2f', r ) );
%         
%         subplot(2,3,5);
%         imshow( -logResponse, [] );
%         hold on;
%         plot(1:size(x,2), hprofile*ones(1,size(x,2)), 'g-');
%         hold off;
%         title( sprintf( 'LoG Response', r ) );
% 
%         subplot(2,3,6);
%         imshow( -lobgResponse, [] );
%         hold on;
%         plot(1:size(x,2), hprofile*ones(1,size(x,2)), 'b-');
%         hold off;
%         title( 'LoBG Response' );
% 
%   Author: Deepak Roy Chittajallu
% 
% 

    p = inputParser;
    p.CaseSensitive( false );
    p.addRequired( 'imInput', @(x) ( isnumeric(x) ) );    
    p.addRequired( 'sigma', @(x) ( isscalar(x) ) ); 
    p.addRequired( 'rho', @(x) ( isscalar(x) ) );
    p.parse( imInput, sigma, rho );    
    
    p.addParamValue( 'spacing', ones( 1, ndims(imInput) ), @(x) (isnumeric(x) && ismember(numel(x), [1,ndims(imInput)])) );
    p.addParamValue( 'borderCondition', 'symmetric', @(x) ( isscalar(x) || (ischar(x) && ismember( lower(x), { 'symmetric', 'replicate', 'antisymmetric' } ) ) ) );
    p.addParamValue( 'UseNormalizedDerivatives', false, @(x) (isscalar(x) && islogical(x)) );
    p.addParamValue( 'UseGPU', false, @(x) (isscalar(x) && islogical(x)) );
    p.parse( imInput, sigma, rho, varargin{:} );

    spacing = p.Results.spacing;
    borderCondition = p.Results.borderCondition;
    flagNormalizeDerivatives = p.Results.UseNormalizedDerivatives;
    flagUseGPU = p.Results.UseGPU;
    
    dims = numel( find(size(imInput) > 1) );

    % adjust sigma according to pixel spacing
    sigmaImsp = sigma * ones( 1, dims ) ./ spacing;

    % Compute the bigaussian kernel
    w = ceil(4 * sigmaImsp);
    xrange = cell(1,dims);
    
    for i = 1:dims
       xrange{i} = -w(i):w(i); 
    end
    
    x = cell(1,dims);
    
    if dims > 1
        [x{:}] = ndgrid( xrange{:} );    
    else
        x{1} = xrange{1};
    end

    rad = zeros(size(x{1}));
    for i = 1:dims
       rad = rad + (x{i} * spacing(i)).^2;
    end    
    rad = sqrt(rad);
    
    normalizeFunc = @(x) ( x ./ sum(x(:)) );
    dimsFunc = @(x) ( numel( find(size(x) > 1) ) );
    guassKernelFunc = @(r,sigma) ( normalizeFunc(exp(-r.^2 / (2*sigma^2))) );
    logKernelFunc = @(r,sigma) ( ((r.^2 - sigma^2)/sigma^4) .* guassKernelFunc(r,sigma) );
    
    fgKernel = logKernelFunc(rad, sigma);
    bgKernel = logKernelFunc(rad+rho*sigma-sigma,rho*sigma);
    
    lobgKernel = fgKernel;
    lobgKernel(rad >= sigma) = rho^2 * bgKernel(rad >= sigma);
    
    if flagNormalizeDerivatives
        lobgKernel = sigma^2 * lobgKernel;
    end
    lobgKernel = lobgKernel - mean( lobgKernel(:) );
    
    % apply in fourier domain
    if dims > 1
        padsize = w;
    else
        padsize = [0, w];
    end
    imPadded = padarrayXT(imInput, padsize, borderCondition);
    imLoBGResponse = convnfft( imPadded, lobgKernel, 'valid', 'UseGPU', flagUseGPU);
    
    % return logKernel if requested
    if nargout > 1
        varargout{1} = lobgKernel;
    end    
end