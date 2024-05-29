function [imSizes] = getImageSizeBatch(filePaths, useParpool)
% Get image sizes without loading the image for a batch of files
% Assume all the images has the same dimensions, 2d or 3d. 
% 
% 
% Author: Xiongtao Ruan (05/22/2024)

if nargin < 2
    useParpool = true;
end

nF = numel(filePaths);

% when there are more than 100 files, use background pool
if useParpool && nF > 100
    p = backgroundPool;
    nworker = p.NumWorkers;

    fs = parallel.FevalFuture;
    for f = 1 : nF
        fs(f) = parfeval(p, @getImageSize, 1, filePaths{f});            
    end
    
    wait(fs, 'finished', nF / nworker * 0.1);
    imSizes = fetchOutputs(fs);
else
    sz = getImageSize(filePaths{1});
    imSizes = zeros(nF, numel(sz));
    imSizes(1, :) = sz;
    for f = 2 : nF
        imSizes(f, :) = getImageSize(filePaths{f});
    end
end

end