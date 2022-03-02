function [region_2, crop_bbox] = crop_subregion_by_intensity(region_2, dimNumThrsh)
% Find the subregion with most rich signal
% 
%  Author: Xiongtao Ruan (02/04/2022)


if nargin < 2 
    dimNumThrsh = 2000;
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
std_2 = std(region_2(:));
region_2_bw = region_2 > mu_2 + 3 * std_2;

for i = 1 : numel(adj_axes)
    ax_i = adj_axes(i);
    % threshold intensity
    fz_2 = squeeze(sum(region_2_bw, setdiff(1 : 3, ax_i)));

    [~, pind] = max(fz_2);

    fa_inds_2 = find(fz_2 > mean(fz_2) + 2 * (0.5 * (std(fz_2) + movstd(fz_2, 11))));

    s = 1;
    t = numel(fa_inds_2);
    sa_2 = 1;
    ta_2 = numel(fz_2);

%     while ta_2 - sa_2 > dimNumThrsh
%         if pind - sa_2 > ta_2 - pind
%             s = s + 1; 
%             sa_2 = fa_inds_2(s);
%         else
%             t = t - 1;
%             ta_2 = fa_inds_2(t);        
%         end        
%     end

    if numel(fz_2) > dimNumThrsh
        if pind - 1 < dimNumThrsh / 2
            sa_2 = 1;
            ta_2 = sa_2 + dimNumThrsh - 1;
        elseif numel(fz_2) - pind < dimNumThrsh / 2
            ta_2 = numel(fz_2);
            sa_2 = ta_2 - dimNumThrsh + 1;
        else
            sa_2 = pind - ceil((dimNumThrsh - 1) / 2) + 1;
            ta_2 = sa_2 + dimNumThrsh - 1;
        end        
    end
    crop_bbox(ax_i) = sa_2;
    crop_bbox(ax_i + 3) = ta_2;
end

region_2 = region_2(crop_bbox(1) : crop_bbox(4), crop_bbox(2) : crop_bbox(5), crop_bbox(3) : crop_bbox(6));


end


