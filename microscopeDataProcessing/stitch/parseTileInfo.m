function [tileInfo, camPatterns] = parseTileInfo(filenames, tilePatterns)
% parse tile patterns for TCXYZ, assume these patterns are separated by
% _ or - or other non-word and non-digits characters, and each pattern has
% a form of word and digit with or without other non-digits characters
% if speficy Cam, Cam should has a pattern of Cam[A-Z], that is, CamA, CamB, ...


nF = numel(filenames);
tileInfo = zeros(nF, 6);

for i = 1 : 5
    tile_pattern_i = tilePatterns{i};
    if isempty(tile_pattern_i)
        continue;
    end
    
    if i == 2
        specifyCam = ~isempty(regexpi(tile_pattern_i, 'Cam[A-Z]', 'once'));
        [tileInds, camInds, camPatterns] = parse_TCXYZ_Info(filenames, tile_pattern_i, specifyCam);
        if specifyCam
            tileInfo(:, 6) = camInds;
        end
    else
        tileInds = parse_TCXYZ_Info(filenames, tile_pattern_i);
    end
    tileInfo(:, i) = tileInds;
end

end


function [tileInds, camInds, camPatterns] = parse_TCXYZ_Info(filenames, tilePattern, specifyCam)
% parse a single pattern for TCXYZ, assume these patterns are separated by
% _ or - or other non-word and non-digits characters, and each pattern has
% a form of word and digit with or without other non-digits characters
% if speficy Cam, Cam should has a pattern of Cam[A-Z], that is, CamA, CamB, ...


if nargin < 3
    specifyCam = false;
end

nF = numel(filenames);
match_ind = regexpi(tilePattern, '\D+', 'once');
if isempty(match_ind)
    error('The tile pattern %s does not contain a non-digit character to specify the pattern!', tilePattern);
end
match_pat = regexpi(tilePattern, '\D+', 'match', 'once');

if match_ind == 1
    pat = sprintf('%s(\\d+)', match_pat);
else
    pat = sprintf('(\\d+)%s', match_pat);
end

tmp = regexpi(filenames, pat, 'tokens');
empty_inds = cellfun(@(x) isempty(x), tmp);

tileInds = zeros(nF, 1);
tileInds(~empty_inds) = cellfun(@(x) str2double(x{1}{1}), tmp(~empty_inds));
tileInds(empty_inds) = NaN;

camInds = [];
camPatterns = [];
if specifyCam
    cpat = 'Cam([A-Z])';
    tmp = regexpi(filenames, cpat, 'tokens');
    empty_inds = cellfun(@(x) isempty(x), tmp);

    camLetters = cell(nF, 1);
    camLetters(~empty_inds) = cellfun(@(x) upper(x{1}{1}), tmp(~empty_inds), 'UniformOutput', false);
    camPatterns = unique(camLetters(~empty_inds));

    camInds = zeros(nF, 1);
    [~, ib] = ismember(camLetters, camPatterns);
    camInds(~empty_inds) = ib;
    camInds(empty_inds) = NaN;
end

end

