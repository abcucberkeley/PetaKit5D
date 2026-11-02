function [im] = unmix_channels_sum(ims, unmixFactors, constBackground, constOffset, dtype, im)
% unmixed image from the channel images: the weighted sum of the channels
% (after subtracting their backgrounds), clipped at 0, plus a constant
% offset, in the given data type
%
% ims: cell array of channel images of the same size
% unmixFactors: unmixing factor of each channel
% constBackground: empty, or the background of each channel (subtracted
%   from the channel before weighting)
% constOffset: empty or a scalar, added after clipping if it is not 0
% dtype: data type of the result, or empty for the partial sum (in single,
%   not clipped and without the offset)
% im (optional): a partial sum (single) to add the channels to, or empty
%   for 0; ims can then be empty. Called as im = unmix_channels_sum(...,
%   im) in a function, MATLAB updates it in place where its code runs.
%
% The mex computes it in one pass (uint8, uint16, single and double images,
% and results); otherwise MATLAB computes it in single. Single channels are
% added to a given partial sum in place by MATLAB, as the original code did
% (the mex returns a new sum); subtract their backgrounds before (in place).


if nargin < 6
    im = [];
end
nc = numel(ims);

supported = @(x) (isa(x, 'uint16') || isa(x, 'uint8') || isa(x, 'single') || isa(x, 'double')) ...
    && isreal(x) && ~issparse(x);
inPlace = nargin > 5 && nc > 0 && all(cellfun(@(x) isa(x, 'single'), ims));
if ~isempty(ims)
    sz = size(ims{1});
else
    sz = size(im);
end
if ~inPlace && (nc > 0 || ~isempty(im)) && all(cellfun(supported, ims)) && all(cellfun(@(x) isequal(size(x), sz), ims)) ...
        && (isempty(dtype) || any(strcmp(dtype, {'uint8', 'uint16', 'single', 'double'}))) ...
        && (isempty(im) || (isa(im, 'single') && isreal(im) && ~issparse(im) && isequal(size(im), sz))) ...
        && isfloat(unmixFactors) && isreal(unmixFactors) && numel(unmixFactors) >= nc && isfloat(constBackground) ...
        && isreal(constBackground) && (isempty(constBackground) || numel(constBackground) >= nc) ...
        && isfloat(constOffset) && isreal(constOffset) && numel(constOffset) <= 1
    bgs = zeros(1, nc);
    if ~isempty(constBackground)
        bgs = double(constBackground(1 : nc));
    end
    offset = 0;
    if ~isempty(constOffset)
        offset = double(constOffset);
    end
    if exist('unmix_channels_sum_mex', 'file') == 3
        try
            im = unmix_channels_sum_mex(ims, double(unmixFactors(1 : nc)), bgs, offset, char(dtype), im);
            return;
        catch ME
            mex_fallback('unmix_channels_sum_mex', ME);
        end
    else
        mex_fallback('unmix_channels_sum_mex');
    end
end

if isempty(im)
    im = zeros(size(ims{1}), 'single');
end
for c = 1 : nc
    im_c = single(ims{c});
    if ~isempty(constBackground)
        im_c = im_c - constBackground(c);
    end
    im = im + im_c * unmixFactors(c);
end
if isempty(dtype)
    return;
end
im = max(0, im);

% add the constant offset after unmixing if it is not 0
if ~isempty(constOffset) && constOffset ~= 0
    im = im + constOffset;
end

im = cast(im, dtype);

end
