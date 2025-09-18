function [fd_ch_groups] = groupImageFilenamesByChannels(filenames, channelPatterns, chInds)
% group images by channels


nF = numel(filenames);
nc = numel(channelPatterns);

% if the third input chInds does not exist, get it from channelPatterns
if nargin < 3
    ch_inds_d = false(nF, nc);
    for c = 1 : nc
        ch_inds_d(:, c) = contains(filenames, channelPatterns{c}, 'IgnoreCase', true) | ...
            contains(filenames, regexpPattern(channelPatterns{c}), 'IgnoreCase', true);
    end
    chInds = arrayfun(@(x) find(ch_inds_d(x, :), 1, 'first'), 1 : size(ch_inds_d, 1));
end

uniq_ch = unique(chInds);
if numel(uniq_ch) ~= nc
    error('Number of channels of the images does not match the number of channel patterns!');
end

fd_ch_groups_cell = cell(nF, 1);
include_flag = false(nF, 1);
group_f_inds = zeros(1, nc);
for f = 1 : nF
    if include_flag(f)
        continue;
    end

    fsn = filenames{f};
    ch_f = chInds(f);
    
    for c = 1 : nc
        if ch_f == c
            group_f_inds(c) = f;
            continue;
        end
        % find the file with largest common substrings from start or end to fsn
        cur_inds = find(~include_flag & chInds(:) == c);
        cur_fsns = filenames(cur_inds);
        common_str_count_mat = zeros(numel(cur_fsns), 1);
        for i = 1 : numel(cur_fsns)
            fsn_i = cur_fsns{i};
            [com_str] = extractCommonStringFromStart(fsn, fsn_i);
            [com_str_1] = extractCommonStringFromStart(flip(fsn), flip(fsn_i));

            common_str_count_mat(i) = numel(com_str) + numel(com_str_1);
        end

        [~, max_ind] = max(common_str_count_mat) ;
        group_f_inds(c) = cur_inds(max_ind);
    end
    include_flag(group_f_inds) = true;
    fd_ch_groups_cell{f} = group_f_inds;
end

fd_ch_groups = cat(1, fd_ch_groups_cell{:});

end


function [com_str] = extractCommonStringFromStart(str_1, str_2)

n1 = numel(str_1);
n2 = numel(str_2);
nc = min(n1, n2);
com_str = str_1(1 : nc);
counter = 0;
for i = 1 : min(n1, n2)
    if str_1(i) == str_2(i)
        com_str(i) = str_1(i);
        counter = i;
    else
        break;
    end
end
com_str = com_str(1 : counter);

end

