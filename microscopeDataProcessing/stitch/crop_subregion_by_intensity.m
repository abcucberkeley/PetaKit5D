function [region_2, crop_bbox] = crop_subregion_by_intensity(region_2, dimNumThrsh, blankNumTrsh)
% Find the subregion with most rich signal
% 
%  Author: Xiongtao Ruan (02/04/2022)
%
%  xruan (05/04/2023): add support to exclude large blank region


if nargin < 2 
    dimNumThrsh = 2000;
end

if nargin < 3
    blankNumTrsh = [];
end

sz_2 = size(region_2);

% if no axis has large length, directly return the region.
crop_bbox = [1, 1, 1, sz_2];
if all(sz_2 <= dimNumThrsh)
    return;
end

[sz_2_adj, adj_axes] = sort(sz_2, 'descend');
adj_axes(sz_2_adj <= dimNumThrsh) = [];

mu_2 = mean(region_2(:));
std_2 = std(single(region_2(:)));
region_2_bw = region_2 > mu_2 + 3 * std_2;

for i = 1 : numel(adj_axes)
    ax_i = adj_axes(i);
    % threshold intensity
    fz_2 = squeeze(sum(region_2_bw, setdiff(1 : 3, ax_i)));

    [~, pind] = max(fz_2);

    s = 1;
    t = numel(fz_2);    

    if ~isempty(blankNumTrsh)
         CC = bwconncomp(fz_2 == 0);
         zero_seq_lengths = cellfun(@numel, CC.PixelIdxList);
         PixelIdxList = CC.PixelIdxList(zero_seq_lengths >= blankNumTrsh);
         PixelIdxList = cat(1, PixelIdxList{:});
         s = max(PixelIdxList(PixelIdxList < pind));
         t = min(PixelIdxList(PixelIdxList > pind));
         if isempty(s)
             s = 1;
         end
         if isempty(t)
             t =  numel(fz_2);
         end
    end

    if t - s + 1 <= dimNumThrsh
        sa_2 = s;
        ta_2 = t;
    else 
        if pind - s + 1 <= dimNumThrsh / 2
            sa_2 = s;
            ta_2 = min(t, sa_2 + dimNumThrsh - 1);
        elseif t - pind + 1 <= dimNumThrsh / 2
            ta_2 = t;
            sa_2 = max(s, ta_2 - dimNumThrsh + 1);
        else
            sa_2 = pind - ceil((dimNumThrsh - 1) / 2) + 1;
            ta_2 = sa_2 + dimNumThrsh - 1;
        end        
    end
    
    crop_bbox(ax_i) = sa_2;
    crop_bbox(ax_i + 3) = ta_2;
end

try 
    region_2 = crop3d_mex(region_2, crop_bbox);
catch ME
    disp(ME);
    region_2 = region_2(crop_bbox(1) : crop_bbox(4), crop_bbox(2) : crop_bbox(5), crop_bbox(3) : crop_bbox(6));
end

end

