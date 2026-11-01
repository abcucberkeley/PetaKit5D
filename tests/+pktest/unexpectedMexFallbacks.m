function [unexpected] = unexpectedMexFallbacks(records, expected, ext)
% the mex fallbacks (records of mex_fallback) that fail a test: of a mex that
% exists ('failed'), or of a missing one that isn't a known gap of the
% platform with mex extension ext (default: this one; pktest.mexPlatformGaps),
% except rejections by design and the mex functions in expected (cell array
% of names the test expects to fall back)

if nargin < 3
    ext = mexext;
end

[gaps, tracked] = pktest.mexPlatformGaps(ext);
keep = false(size(records));
for i = 1 : numel(records)
    r = records(i);
    knownGap = strcmp(r.kind, 'missing') && (~tracked || any(strcmp(r.mex, gaps)));
    keep(i) = ~strcmp(r.kind, 'designed') && ~knownGap && ~any(strcmp(r.mex, expected));
end
unexpected = records(keep);

end
