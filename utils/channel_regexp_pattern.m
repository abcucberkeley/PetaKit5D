function [pat] = channel_regexp_pattern(channelPattern)
% the regular expression of a channel pattern for contains(), as
% regexpPattern(channelPattern) gives it. A pattern that is not a valid
% regular expression (e.g., with an unmatched '(', '[' or '{', as 'Cam(A')
% is matched as text instead, as the literal test next to it does, instead
% of failing.

try
    pat = regexpPattern(channelPattern);
catch
    pat = channelPattern;
end

end
