function [imout] = imresize3_average(im, resampleFactor)
% resize by averaging the data

rs = resampleFactor;
sz = size(im, [1, 2, 3]);

if any(rem(sz, rs) ~= 0)
    pad_size = ceil(sz ./ rs) .* rs - sz;
    im = padarray(im, pad_size, 0, 'post');
    sz_1 = size(im);
    im_c = false(sz_1);
    im_c(1 : sz(1), 1 : sz(2), 1 : sz(3)) = true;
    
    im = single(im);
    im = reshape(im, rs(1), sz_1(1) / rs(1), rs(2), sz_1(2) / rs(2), rs(3), sz_1(3) / rs(3)); 
    im_c = reshape(im_c, rs(1), sz_1(1) / rs(1), rs(2), sz_1(2) / rs(2), rs(3), sz_1(3) / rs(3)); 
    
    imout = squeeze(sum(im, [1, 3, 5]));
    im_c = squeeze(sum(im_c, [1, 3, 5]));
    imout = imout ./ im_c;
else
    im = single(im);
    im = reshape(im, rs(1), sz(1) / rs(1), rs(2), sz(2) / rs(2), rs(3), sz(3) / rs(3));
    imout = squeeze(mean(im, [1, 3, 5]));
end

end

