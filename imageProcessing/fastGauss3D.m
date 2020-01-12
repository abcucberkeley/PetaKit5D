function out=fastGauss3D(img,sigma,fSze,correctBorder,filterMask,reduceNanImage)
% fastGauss3D	apply a 2/3 dimensional gauss filter
%
%    SYNOPSIS out=fastGauss3D(img,sigma,fSze,correctBorder,filterMask)
%
%    INPUT:   img      2 or 3-dimensional data
%             sigma    of gauss filter. Supply single sigma or sigma for
%                      every dimension
%             fSze     (optional) size of the gauss mask [sizeX sizeY sizeZ]
%                          (odd size required for symmetric mask!)
%                       If empty, odd mask with +/- 4*sigma is used.
%             correctBorder (optional) if 1, border effects of the filtering are
%                              lessened. If 2, old version of correction is
%                              used. Default: 1. Note: correctBorder = 2
%                              requires 3D image.
%             filterMask    (optional) supply your own mask. In this case,
%                              sigma and fSze don't have to be supplied. If
%                              the filter mask is a cell, the data will be
%                              sequentially filtered by the contents of the
%                              cell. In that case, fSze is required.
%             reduceNanImage (optional) if 1, nan-containing areas are
%                              clipped before filtering to increase speed.
%                              Default: 0
%
%    OUTPUT:  out      filtered data
%
% c: 13/03/01 dT
% revamped by jonas

%===============
%% check input
%===============
if nargin < 1 || isempty(img);
    error('Please pass nonempty image to fastGauss3D')
end
% set default correctBorder
if nargin < 4 || isempty(correctBorder)
    correctBorder = 1;
end

% read dimensionality
dims = sum(size(img)>0);

% check for filterMask
if nargin < 5 || isempty(filterMask)
    % in this case, we need fSze etc
    if nargin < 2 || isempty(sigma) || ~any(length(sigma)==[1 dims])
        error('please supply nonempty sigma of correct dimensionality')
    end
    if length(sigma) == 1
        sigma = sigma * ones(1,dims);
    end
    % check for filterSize
    if nargin < 3 || isempty(fSze)
        fSze = roundOddOrEven(sigma(1:dims)*4,'odd','inf');
    end
    % create gaussMask
    switch dims
        case 2
            filterMask=GaussMask2D(sigma,fSze,[],1,[],[],1);
        case 3
            filterMask=GaussMask3D(sigma,fSze,[],1,[],[],1);
    end
else
    if isempty(fSze)
        if iscell(filterMask)
            error('if you supply a separated filter, you will need to supply the filter size!')
        else
            fSze = size(filterMask);
        end
    end
end

if nargin < 6 || isempty(reduceNanImage)
    reduceNanImage = false;
end


% add border to image
convnOpt = 'same';
nanMask = []; % mask indicating all the NaNs in the original image;
if correctBorder == 2
    if dims < 3
        warning('FASTGAUSS3D:WRONGDIMENSION',...
            'Cannot use old correctBorder if not 3D image. Using new correctBorder instead')
        correctBorder = 1;
    else
        addBorderOld;
        convnOpt = 'valid';
    end
end
if correctBorder == 1
    %correct for border effects. Create nanMask first. It's a logical mask
    %that takes up 1/8th the space of the original. As long as more than
    %1/8th of the image contains Nan, it's smaller than a list of indices
    nanMask = isnan(img); 
    fullMask = [];
    if reduceNanImage && any(nanMask(:))
        goodRCZcell = cell(dims,1);
        nan2 = all(nanMask,3);
        goodRCZcell{1} = ~all(nan2,2);
        goodRCZcell{2} = ~all(nan2,1);
        clear nan2
        if dims > 2
            goodRCZcell{3} = ~all(all(nanMask,1),2);
        end
        nanRatio = 1-prod(cellfun(@(x)(sum(x)),goodRCZcell))/numel(nanMask);
        if nanRatio > 0.05 %-- from a little bit of testing it looks like you get about nanRatio*0.4 reduction in time
            img = img(goodRCZcell{:});
            fullMask = nanMask;
            nanMask = nanMask(goodRCZcell{:});
        end
    end
        
    % pass nanMask so that we don't need to recalc again
    img = addBorder(img,floor(fSze/2),nanMask);
    convnOpt = 'valid';
end

switch dims
    case 2
        % Convolve matrices
        if iscell(filterMask)
            for i=1:length(filterMask)
                if length(filterMask{i}) > 1
                    img = conv2(img,filterMask{i},convnOpt);
                end
            end
            out = img;
        else
            out=conv2(img,filterMask,convnOpt);
        end
    case 3
        % Convolve matrices
        if iscell(filterMask)
            for i=1:length(filterMask)
                if length(filterMask{i}) > 1
                    img = convn(img,filterMask{i},convnOpt);
                end
            end
            out = img;
        else
            out=convn(img,filterMask,convnOpt);
        end
end


if ~isempty(nanMask)
    if reduceNanImage && ~isempty(fullMask)
        % rebuild full image
        img = out;
        out = NaN(size(fullMask));
        out(~fullMask) = img(~nanMask);
    else
        % simply ensure that NaNs are in the right place
    out(nanMask) = NaN;
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function addBorderOld
        %ADDBORDER adds a border around 3D images for filtering
        %adds a border of halFsze around the img, and fill it with pixels whose
        %value is computed as follows:
        %for all 3 dimensions i (=6 sides) of the array, the nanmedian of the first or last
        %halFsze(i)+1 pixels in this dimension, respectively, is computed. The
        %edges and corners are filled with means of the adjacent sides or edges,
        %respectively
        %c: 030514 jonas

        %init
        halFsze = zeros(1,3);
        halFsze(1:length(fSze)) = floor(fSze/2); %half filter size
        oldSize = size(img);
        if length(oldSize) < 3
            error('image too small for 3D filtering: no 3rd dimension!')
        elseif any(oldSize < halFsze+1)
            error('image too small for 3D filtering: filter is larger than image!')
        end

        newImg = zeros(oldSize+2*halFsze);
        newImg(halFsze(1)+1:end-halFsze(1),halFsze(2)+1:end-halFsze(2),halFsze(3)+1:end-halFsze(3)) = img;

        %alternatively, this could be made with feval and for-loop, but it would be a
        %chore, too

        %2 sides first direction
        inSide11 = img(1:halFsze(1)+1,:,:);
        side11val = nanmedian(inSide11(:));
        newImg(1:halFsze(1),halFsze(2)+1:end-halFsze(2),halFsze(3)+1:end-halFsze(3)) = side11val;

        inSide12 = img(end-halFsze(1):end,:,:);
        side12val = nanmedian(inSide12(:));
        newImg(end-halFsze(1)+1:end,halFsze(2)+1:end-halFsze(2),halFsze(3)+1:end-halFsze(3)) = side12val;
        %2 sides second direction
        inSide21 = img(:,1:halFsze(2)+1,:);
        side21val = nanmedian(inSide21(:));
        newImg(halFsze(1)+1:end-halFsze(1),1:halFsze(2),halFsze(3)+1:end-halFsze(3)) = side21val;

        inSide22 = img(:,end-halFsze(2):end,:);
        side22val = nanmedian(inSide22(:));
        newImg(halFsze(1)+1:end-halFsze(1),end-halFsze(2)+1:end,halFsze(3)+1:end-halFsze(3)) = side22val;

        %2 sides third direction
        inSide31 = img(:,:,1:halFsze(3)+1);
        side31val = nanmedian(inSide31(:));
        newImg(halFsze(1)+1:end-halFsze(1),halFsze(2)+1:end-halFsze(2),1:halFsze(3)) = side31val;

        inSide32 = img(:,:,end-halFsze(3):end);
        side32val = nanmedian(inSide32(:));
        newImg(halFsze(1)+1:end-halFsze(1),halFsze(2)+1:end-halFsze(2),end-halFsze(3)+1:end) = side32val;

        %4 edges parallel to first direction
        edge2131val = nanmean([side21val,side31val]);
        newImg(:,1:halFsze(2),1:halFsze(3)) = edge2131val;

        edge2231val = nanmean([side22val,side31val]);
        newImg(:,end-halFsze(2)+1:end,1:halFsze(3)) = edge2231val;

        edge2132val = nanmean([side21val,side32val]);
        newImg(:,1:halFsze(2),end-halFsze(3)+1:end) = edge2132val;

        edge2232val = nanmean([side22val,side32val]);
        newImg(:,end-halFsze(2)+1:end,end-halFsze(3)+1:end) = edge2232val;
        % 4 edges parallel to second direction
        edge1131val = nanmean([side11val,side31val]);
        newImg(1:halFsze(1),:,1:halFsze(3)) = edge1131val;
        edge1231val = nanmean([side12val,side31val]);
        newImg(end-halFsze(1)+1:end,:,1:halFsze(3)) = edge1231val;
        edge1132val = nanmean([side11val,side32val]);
        newImg(1:halFsze(1),:,end-halFsze(3)+1:end) = edge1132val;
        edge1232val = nanmean([side12val,side32val]);
        newImg(end-halFsze(1)+1:end,:,end-halFsze(3)+1:end) = edge1232val;
        % 4 edges parallel to third direction
        edge1121val = nanmean([side11val,side21val]);
        newImg(1:halFsze(1),1:halFsze(2),:) = edge1121val;

        edge1221val = nanmean([side12val,side21val]);
        newImg(end-halFsze(1)+1:end,1:halFsze(2),:) = edge1221val;

        edge1122val = nanmean([side11val,side22val]);
        newImg(1:halFsze(1),end-halFsze(2)+1:end,:) = edge1122val;
        edge1222val = nanmean([side12val,side22val]);
        newImg(end-halFsze(1)+1:end,end-halFsze(2)+1:end,:) = edge1222val;
        %corner 000
        newImg(1:halFsze(1),1:halFsze(2),1:halFsze(3)) = nanmean([edge2131val,edge1131val,edge1121val]);
        %corner 100
        newImg(end-halFsze(1)+1:end,1:halFsze(2),1:halFsze(3)) = nanmean([edge2131val,edge1231val,edge1221val]);
        %corner 010
        newImg(1:halFsze(1),end-halFsze(2)+1:end,1:halFsze(3)) = nanmean([edge2231val,edge1131val,edge1122val]);

        %corner 110
        newImg(end-halFsze(1)+1:end,end-halFsze(2)+1:end,1:halFsze(3)) = nanmean([edge2231val,edge1231val,edge1222val]);

        %corner 001
        newImg(1:halFsze(1),1:halFsze(2),end-halFsze(3)+1:end) = nanmean([edge2132val,edge1132val,edge1121val]);

        %corner 101
        newImg(end-halFsze(1)+1:end,1:halFsze(2),end-halFsze(3)+1:end) = nanmean([edge2132val,edge1232val,edge1221val]);

        %corner 011
        newImg(1:halFsze(1),end-halFsze(2)+1:end,end-halFsze(3)+1:end) = nanmean([edge2232val,edge1132val,edge1122val]);
        %corner 111
        newImg(end-halFsze(1)+1:end,end-halFsze(2)+1:end,end-halFsze(3)+1:end) = nanmean([edge2232val,edge1232val,edge1222val]);


        % assign img
        img = newImg;
    end
end
