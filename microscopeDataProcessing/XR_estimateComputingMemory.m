function [estRequiredMemory, estRequiredGPUMem, rawImageSize, imSize] = XR_estimateComputingMemory(filePath, varargin)
% estimate memory requirement for the computing
% The unit is Gb
% 
% Author: Xiongtao Ruan (02/27/2020)
% xruan (02/27/2020): fix issue for tall image for deskew. 
% xruan (08/28/2020): add support of input of image size 


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('filePath'); 
ip.addOptional('steps', {'deskew', 'rotate', 'deconvolution'}); 
ip.addParameter('imSize', [], @(x) isnumeric(x) && (isempty(x) || numel(x) == 3)); 
ip.addParameter('memFactors', [10, 5, 10]); 
ip.addParameter('cudaDecon', ~false, @islogical);
ip.addParameter('GPUMemFactor', 1.5); 
ip.addParameter('GPUMaxMem', 12, @isnumeric);

ip.parse(filePath, varargin{:});

warning('off', 'all');

steps = ip.Results.steps;
imSize = ip.Results.imSize;
memFactors = ip.Results.memFactors;
cudaDecon = ip.Results.cudaDecon;
GPUMemFactor = ip.Results.GPUMemFactor;

if ~exist(filePath, 'file')
    if isempty(imSize)
        error('File %s does not exist!', filePath);        
    end
else
    imSize = getImageSize(filePath);    
end

% in double size
rawImageSize = prod(imSize) * 8 / 1024^3;

estRequiredMemory = zeros(numel(steps), 1);

if contains(steps, 'deskew', 'IgnoreCase', true)
    ind = strcmpi(steps, 'deskew');
    % use 300 slices as threshold
    estRequiredMemory(ind) = rawImageSize * memFactors(1) * max(1, (imSize(3) / 300) ^ 2);
end

if contains(steps, 'rotate', 'IgnoreCase', true)
    ind = strcmpi(steps, 'rotate');
    estRequiredMemory(ind) = rawImageSize * memFactors(2);
end

if contains(steps, 'deconvolution', 'IgnoreCase', true)
    ind = strcmpi(steps, 'deconvolution');
    estRequiredMemory(ind) = rawImageSize * memFactors(3);
    estRequiredGPUMem = NaN;
    if cudaDecon
        estRequiredGPUMem = rawImageSize * GPUMemFactor;
    end
end

end