function [idx, sigma_r_thrsh, sure_level] = filter_bad_points_in_detection(pM, sigma_r, ResidualSigmaThresh, BackgroundInfo, vol, A, c)
% identify bad points in the detection based on multiple criterion:
% 1. isolated points (second nearest distace)
% 2. multiple peaks in the left side in A, c, sigma_r
% 
% Author: Xiongtao Ruan 12/06/2019


background_sigma = BackgroundInfo(3);
ResidualSigmaThresh_device = BackgroundInfo(4);

% define if the background sigma is too large, in which case the image is
% noisy
large_bg_std = false;
if ResidualSigmaThresh_device < 5
    if ResidualSigmaThresh_device * 2 < background_sigma % SCMOS
        large_bg_std = true;
    end
else
    if ResidualSigmaThresh_device * 3 < background_sigma % EMCCD
        large_bg_std = true;
    end        
end

sure_level = 0;

if numel(sigma_r) < 3
    idx = true;
    sigma_r_thrsh = -1;
    return;
end

% identify background points
sigma_high_low_ratio = prctile(sigma_r, 95) / prctile(sigma_r, 5);

sigma_ratio_thresh = 2.0;
sigma_r_thrsh = -1;

if sigma_high_low_ratio > sigma_ratio_thresh
    if ResidualSigmaThresh_device < 5
        interval_bound = 4; % SCMOS
    else
        interval_bound = 15; % EMCCD
    end
    % in case of there are very large stds, bound the maximum to 20 fold of
    % std device or 95% quantile
    if prctile(sigma_r, 99) > 20 * ResidualSigmaThresh_device
        sigma_r_ub = max(20 * ResidualSigmaThresh_device, prctile(sigma_r, 95));
    else
        sigma_r_ub = prctile(sigma_r, 99);
    end
    
    bin_intervel = min(max(std(sigma_r) / 10, (prctile(sigma_r, 95) - prctile(sigma_r, 5)) / 90), interval_bound);
    bin_num_0 = ceil((sigma_r_ub - min(sigma_r)) / bin_intervel);
    % Freedman-Diaconis
    bin_num_1 = ceil((sigma_r_ub - min(sigma_r)) / (2 * (prctile(sigma_r, 75) - prctile(sigma_r, 25)) / numel(sigma_r) ^ (1 / 3))); 
    bin_num = min(ceil(numel(sigma_r) / 5), max([60, bin_num_0, bin_num_1]));
    [N, edges] = histcounts(sigma_r(sigma_r < sigma_r_ub), bin_num);
    
    if numel(N) > 2
        [pkds, locs] = findpeaks(N, 'minpeakprominence', 5);
        % in case there are lots of local maximum, reduce local maximum by
        % reducing bin number. 
        if numel(locs) > 20 || sum(diff(locs(locs < bin_num / 4)) < 5) > 3 || (sum(N(diff(N) == 0) >  prctile(N, 25)) > 2 && numel(sigma_r) > 100)
            bin_num = ceil(bin_num * 0.9);
            [N, edges] = histcounts(sigma_r(sigma_r < sigma_r_ub), bin_num);
            [pkds, locs] = findpeaks(N, 'minpeakprominence', 5);

            while numel(locs) > 20 || sum(diff(locs(locs < bin_num / 4)) < 5) > 3 || (sum(N(diff(N) == 0) >  prctile(N, 25)) > 2 && numel(sigma_r) > 100)
                bin_num = ceil(bin_num * 0.9);
                [N, edges] = histcounts(sigma_r(sigma_r < sigma_r_ub), bin_num);
                [pkds, locs] = findpeaks(N, 'minpeakprominence', 5);
            end
            
            [pkds, locs] = findpeaks(N, 'minpeakprominence', 5);
        end
        
        edge_c = (edges(1 : end - 1) + edges(2 : end)) / 2;
        
        % in case some local maximum in the lower side is not very prominent
        if ~isempty(locs)
            % if ~any(edge_c(locs) < prctile(sigma_r, 25))
            [~, max_low_ind] = max(N(edge_c < prctile(sigma_r, 30)));
            if ~any(ismember(max(1, max_low_ind - 1) : max_low_ind + 1, locs))
                locs = [max_low_ind, locs];
            end
            
            if N(1) > N(2) && all(locs ~= 1)
                locs = [1, locs];
            end
        end
        
        % if there are a narrow gap in the lower side, smooth the histogram via median filter
        if numel(locs) > 3 && bin_num / 5 > locs(4)
            loc_lb = locs(4);
        else
            loc_lb = bin_num / 5;
        end
        if sum(locs < loc_lb) > 1 && any(diff(locs(locs < loc_lb)) < 4) 
            N_med = medfilt1(N);
            N_diff = N_med - N;
            % if there is dramatic change in local maximum, use moving
            % average filter instead
            if any(abs(N_diff(locs < loc_lb) ./ N(locs < loc_lb)) > 0.25 & N(locs < loc_lb) > median(N(locs)))
                B = [0.1, 0.8, 0.1];
                N_ma = conv(N, B, 'same');
                N = N_ma;
            else
                N = N_med;
            end
            [pkds, locs] = findpeaks(N, 'minpeakprominence', 5);                
        end
    else
        idx = true(1, numel(sigma_r));
        return;
    end
    
    % in case some local maximum in the lower side is not very prominent
    if ~isempty(locs)
        % if ~any(edge_c(locs) < prctile(sigma_r, 25))
        % iterative add low side not prominent local maximum
        for i = 1 : 3
            [~, max_low_ind] = max(N(edge_c < prctile(sigma_r, 40 - i * 10)));
            if ~isempty(max_low_ind) && ~any(ismember(max(1, max_low_ind - 1) : max_low_ind + 1, locs)) ...
                    && (max_low_ind == 1 || (N(max_low_ind) >= N(max_low_ind - 1) && N(max_low_ind) >= N(max_low_ind + 1)))
                locs = [max_low_ind, locs];
            end
        end
        locs = sort(locs);
        
        if N(1) > N(2) && all(locs ~= 1)
            locs = [1, locs];
        end
    end    
    
    if numel(locs) > 1
        % the first peak need to be in the very left side, and contains lots of
        % points, and contains at least three bins (in case of the left side
        % local maximum bin contains lots of points)
        left_loc_inds = edge_c(locs) <= min(max(prctile(sigma_r, 30), edge_c(4)), edge_c(round(bin_num / 5)));
        % high_peak_inds = N(locs) > 0.01 * numel(sigma_r);
        if numel(locs) > 2
            high_peak_inds = N(locs) >= min(prctile(N, 85), max(prctile(N, 70), median(N(locs))));
            % if the first loc is tall enough, set it as true
            if ~high_peak_inds(1) && N(locs(1)) / sum(N) > 0.01
                 high_peak_inds(1) = true;
            end                
        else
            high_peak_inds = true(1, 2);
        end
        % bound on the threshold.
        if large_bg_std
            bounded_sigma_inds = edge_c(locs) < 2.5 * background_sigma;
        else
            bounded_sigma_inds = edge_c(locs) < 2 * max(background_sigma, ResidualSigmaThresh_device);
        end

        left_high_ind = find(left_loc_inds & high_peak_inds & bounded_sigma_inds);

        if any(left_high_ind)
            left_high_ind_all = left_high_ind;
%             if numel(locs) < 3               
%                 left_high_ind = left_high_ind(1);
%             else
%                 if left_high_ind(end) == numel(locs)
%                     left_high_ind = left_high_ind(end - 1);
%                 else
%                     left_high_ind = left_high_ind(end);
%                 end
%             end
            if left_high_ind(end) == numel(locs)
                if numel(left_high_ind) > 1
                    left_high_ind = left_high_ind(end - 1);
                else
                    left_high_ind = left_high_ind - 1;
                end
            else
                if numel(left_high_ind) > 1 
                    [~, max_ind] = max(N(locs(left_high_ind)));
                    if N(locs(left_high_ind(max_ind))) / N(locs(left_high_ind(end))) > 1.5
                        left_high_ind = left_high_ind(max_ind);
                    else
                        left_high_ind = left_high_ind(end);
                    end
                else
                    left_high_ind = left_high_ind(end);                    
                end
            end                  
            
            chosen_locs = locs([left_high_ind, left_high_ind + 1]);
            [~, min_ind] = min(N(chosen_locs(1) : chosen_locs(2)));
            min_ind = min_ind + chosen_locs(1) - 1;
            
            % in case the minumum is because of noise, check whether it is
            % a deep valley, if not, do not use it. 
            if N(min_ind) / N(locs(left_high_ind)) < 0.6
                sigma_r_thrsh = edge_c(min_ind);
                if N(min_ind) / N(locs(left_high_ind)) < 0.3
                    sure_level = 1;
                end
                if N(min_ind) / N(locs(left_high_ind)) < 0.15 && ...
                    ((ResidualSigmaThresh_device < 5 && sigma_r_thrsh < 2 * ResidualSigmaThresh_device) || ...  % SCMOS
                    (ResidualSigmaThresh_device > 5 && sigma_r_thrsh < 3 * ResidualSigmaThresh_device)) && ...  % EMCCD
                    edge_c(locs(left_high_ind)) < 1.5 * ResidualSigmaThresh_device
                    sure_level = 3;
                end
                
                % if the local maximum is very tall and sigma_r_thrsh is very large
                % but sure level is not 3, we can reduce the threshold by 
                % using the boundary of this local maximum. 
                
                if sure_level ~= 3 
                    % in case there are multiple local maximum, use the
                    % tallest one. 
                    if numel(left_high_ind_all) > 1 && max(N(locs(left_high_ind_all))) ~= N(locs(left_high_ind))
                        [~, max_ind] = max(N(locs(left_high_ind_all)));
                        h_loc = locs(left_high_ind_all(max_ind));
                    else
                        h_loc = locs(left_high_ind);
                    end
                    if N(h_loc) / sum(N) > 0.08 || sum(N(max(1, h_loc -1) : min(bin_num, max(h_loc + 1, 3)))) / sum(N) > 0.12 ...
                            || sum(maxk(N(max(1, h_loc -1) : min(bin_num, max(h_loc + 1, 3))), 2)) / sum(N) > 0.09
                        loc_half = [];
                        n_factor = 4;
                        n_factor_lb = 2;
                        while isempty(loc_half) && n_factor >= n_factor_lb
                            loc_half = find(N < N(h_loc) / n_factor & (1 : bin_num) > h_loc & (1 : bin_num) < locs(left_high_ind + 1), 1, 'first');
                            n_factor = n_factor - 0.5;
                        end
                            
                        if ~isempty(loc_half)    
                            h_ind = max(loc_half, h_loc + 1);
                        else
                            h_ind = h_loc + 1;
                        end
                        sigma_r_thrsh_alt = edges(h_ind + 1);
                        
                        if ((ResidualSigmaThresh_device < 5 && sigma_r_thrsh_alt < 2 * ResidualSigmaThresh_device) || ...  % SCMOS
                            (ResidualSigmaThresh_device > 5 && sigma_r_thrsh_alt < 3 * ResidualSigmaThresh_device)) && ...  % EMCCD
                            edge_c(h_loc) < 1.5 * ResidualSigmaThresh_device && N(h_loc) / N(h_ind) > 2.5
                            sigma_r_thrsh = sigma_r_thrsh_alt;
                            sure_level = 3;
                        end
                        % if the peak is very tall, force std thrsh to the
                        % alternative one. 
                        if sure_level ~= 3 &&  ...
                                (edge_c(h_loc) < max(1.5 * ResidualSigmaThresh_device, background_sigma) || N(h_loc) / sum(N) > 0.1) ...
                                && N(h_loc) / N(h_ind) > 5 && ...
                                ((ResidualSigmaThresh_device < 5 && sigma_r_thrsh_alt < 5 * max(ResidualSigmaThresh_device, background_sigma)) || ...  % SCMOS
                                (ResidualSigmaThresh_device > 5 && sigma_r_thrsh_alt < 10 * max(ResidualSigmaThresh_device, background_sigma)))  % EMCCD
                            sigma_r_thrsh = sigma_r_thrsh_alt;
                            sure_level = 3;
                        end
                        % if the peak is tall enought, force to increase
                        % the sigma_r_thrsh
                        if sure_level ~= 3 && ...
                                (edge_c(h_loc) < max(1.5 * ResidualSigmaThresh_device, background_sigma) || N(h_loc) / sum(N) > 0.1) ...
                                && N(h_loc) / N(h_ind) > 3 && ...
                                ((ResidualSigmaThresh_device < 5 && sigma_r_thrsh_alt < 5 * max(ResidualSigmaThresh_device, background_sigma)) || ...  % SCMOS
                                (ResidualSigmaThresh_device > 5 && sigma_r_thrsh_alt < 10 * max(ResidualSigmaThresh_device, background_sigma)))  % EMCCD                                
                            if ResidualSigmaThresh_device < 5
                                sigma_r_thrsh = min(1.5 * ResidualSigmaThresh_device, background_sigma);
                            else
                                sigma_r_thrsh = min(2 * ResidualSigmaThresh_device, background_sigma);
                            end
                            sure_level = 2;
                        end
                    end
                end
            end
        end
    end
end

% estimate isolated scatter points which are probably background points. 
% [local_density] = estimate_local_properties(pM);

[~, dist_mat] = knnsearch(pM, pM, 'k', 3);
% r_radius = median(dist_mat(:, 3)) * 2;
r_radius = mean(dist_mat(:, 3)) + 3 * std(dist_mat(:, 3));
[idx_cell] = rangesearch(pM, pM, r_radius);
local_density = cellfun(@numel, idx_cell);
local_density = local_density';

rm_iso_pt = local_density == 1;

rm_sigma_r = sigma_r < sigma_r_thrsh;

if sum(rm_sigma_r) > 0
    if sum(rm_iso_pt) / numel(rm_iso_pt) > 0.05 && sum(rm_iso_pt .* rm_sigma_r) / sum(rm_iso_pt) > 0.8 && sure_level < 2
        sure_level = 2;
    end
else
    if sum(rm_iso_pt) / numel(rm_iso_pt) > 0.02
        p_i = 2.5;
        sigma_r_thrsh_1 = prctile(sigma_r(rm_iso_pt), p_i);
        rm_sigma_r_1 = sigma_r < sigma_r_thrsh_1;
        while sum(rm_iso_pt .* rm_sigma_r_1) / min(sum(rm_iso_pt), sum(rm_sigma_r_1)) > 0.9
            sigma_r_thrsh_1 = prctile(sigma_r(rm_iso_pt), p_i);
            rm_sigma_r_1 = sigma_r < sigma_r_thrsh_1;
            p_i = p_i + 2.5;
            if p_i >= 100
                break;
            end
        end
        sigma_r_thrsh = sigma_r_thrsh_1;
        rm_sigma_r = rm_sigma_r_1;
        if sum(rm_iso_pt .* rm_sigma_r) / sum(rm_iso_pt) > 0.85 && sure_level < 2
            sure_level = 2;
        end
    end   
end


if sure_level >= 2 || (sure_level == 1 && sigma_r_thrsh > ResidualSigmaThresh && ~large_bg_std && sigma_r_thrsh ...
        < max(2 * max(background_sigma, ResidualSigmaThresh_device), 2.5 * min(background_sigma, ResidualSigmaThresh_device)))
    idx = ~(rm_sigma_r);
elseif sure_level < 2
    idx = sigma_r > ResidualSigmaThresh;
end

% keep points with high local density (with high enought sigma_r) and remove isolated points
% idx(~idx) = local_density(~idx) > 3;
% idx(idx) = local_density(idx) > 1;

local_density_thresh_0 = 8;
if sure_level < 3
    local_density_thresh = 3;
    if sigma_r_thrsh == -1
        local_density_thresh = 2;
    end
else
    local_density_thresh = 4;
end    

if sure_level >= 2
    idx_0 = (idx | local_density > local_density_thresh_0 | (local_density > local_density_thresh & sigma_r > ...
        max([0.8 * sigma_r_thrsh, 0.9 * background_sigma, min(background_sigma, ResidualSigmaThresh_device)]))) & (local_density > 1);
else
    if 0.8 * sigma_r_thrsh * 1.5 < max(0.9 * background_sigma, min(background_sigma, ResidualSigmaThresh_device))
        idx_0 = (idx | local_density > local_density_thresh_0 | (local_density > local_density_thresh & sigma_r > ...
            max([0.8 * sigma_r_thrsh, 0.9 * background_sigma, min(background_sigma, ResidualSigmaThresh_device)]))) & (local_density > 1);
    else
        idx_0 = (idx | local_density > local_density_thresh_0 | (local_density > local_density_thresh & sigma_r > ...
            max([0.9 * background_sigma, min(background_sigma, ResidualSigmaThresh_device)]))) & (local_density > 1);        
    end
end

% further check if we keep too many noise point because of dense neighbors.

if sure_level >= 2
    if sum(~idx_0) / sum(sigma_r < sigma_r_thrsh) < 0.5 && sum(sigma_r < sigma_r_thrsh) > 100
        local_density_thresh_0 = max(9, min(10, median(local_density(sigma_r < sigma_r_thrsh))));
    end
    idx = (idx | (local_density > local_density_thresh_0) ...
        | (local_density > local_density_thresh & sigma_r > ...
        max([0.8 * sigma_r_thrsh, 0.9 * background_sigma, min(background_sigma, ResidualSigmaThresh_device)]))) & (local_density > 1);        
else
    if sure_level == 1 && sum(~idx_0) / sum(sigma_r < sigma_r_thrsh) < 0.5 && sum(sigma_r < sigma_r_thrsh) > 100
        local_density_thresh_0 = max(9, min(10, median(local_density(sigma_r < sigma_r_thrsh))));
    end    
    if 0.8 * sigma_r_thrsh * 1.5 < max(0.9 * background_sigma, min(background_sigma, ResidualSigmaThresh_device))
        idx = (idx | (local_density > local_density_thresh_0) ...
            | (local_density > local_density_thresh & sigma_r > ...
            max([0.8 * sigma_r_thrsh, 0.9 * background_sigma, min(background_sigma, ResidualSigmaThresh_device)]))) & (local_density > 1);
    else
        idx = (idx | (local_density > local_density_thresh_0) ...
            | (local_density > local_density_thresh & sigma_r > ...
            max([0.9 * background_sigma, min(background_sigma, ResidualSigmaThresh_device)]))) & (local_density > 1);        
    end
end


end


