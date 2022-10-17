function [img] = XR_locmax3d(img, wdims, ClearBorder)
% local maximum based on image dilation 
% The function is copied from Francis Agust locamax3d.m
% The idea is from https://stackoverflow.com/questions/1856197/how-can-i-find-local-maxima-in-an-image-in-matlab?rq=1
% Here we treat NaN values as -Inf. 
% 
% Author: Xiongtao Ruan 10/12/2019


% ip = inputParser;
% ip.CaseSensitive = false;
% ip.addRequired('img', @isnumeric);
% ip.addRequired('wdims', @isnumeric);
% ip.addParameter('ClearBorder', true, @islogical);
% ip.parse(img, wdims, varargin{:});

if nargin < 3
    ClearBorder = false;
end    

if numel(wdims)==1
    wx = wdims;
    wy = wdims;
    wz = wdims;
    if mod(wx,2)==0 || mod(wy,2)==0 || mod(wz,0)==0
        error('Mask dimensions must be odd integers');
    end
elseif numel(wdims)==3
    wx = wdims(1);
    wy = wdims(2);
    wz = wdims(3);
    if mod(wx,2)==0 || mod(wy,2)==0 || mod(wz,0)==0
        error('Mask dimensions must be odd integers');
    end    
end

hwx = (wx - 1) / 2;
hwy = (wy - 1) / 2;
hwz = (wz - 1) / 2;

% define kernel 
se = true(wx, wy, wz);
se(hwy + 1, hwx + 1, hwz + 1) = false;

% img_1 = img;
is_nan_mat = isnan(img);
% if any(is_nan_mat, 'all')
%      img_1(is_nan_mat) = -1e8;
% end

if any(is_nan_mat, 'all')
     img = replace_nan_with_value(img, -1e8);
     % img_1(is_nan_mat) = -1e8;
end    

% img_3 = replace_nan_with_value(img, -1e8);

% 11/13/2019 xruan use builtin function for imdilate
B = imdilate(img, se);
% minmax = [-inf; inf];
% morphFunc = 'imdilate';
% B = builtin('_morphmex_halide', img_1, minmax, se,  morphFunc);
% bw = img_1 > B;

% if max. filter response is equal to input, point is a local maximum
img = img .* (img > B); 

% set borders to zero
if ClearBorder
    b = hwx;
    img(:,[1:b, end-b+1:end],:) = 0;
    b = hwy;
    img([1:b, end-b+1:end],:,:) = 0;
    b = hwz;
    img(:,:,[1:b, end-b+1:end]) = 0;
end

end

