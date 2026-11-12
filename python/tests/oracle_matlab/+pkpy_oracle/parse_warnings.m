function [w] = parse_warnings(txt, lastMsg, lastId)
% the warnings in captured command window text (verbose on, backtrace off):
% a cell array of structs with identifier and message. Same rules as
% petakit5d.testing.oracle.parse_matlab_warnings.
%
% MATLAB writes a warning as [<backspace>Warning: message]<backspace>, with
% the line (Type "warning off <id>" to suppress this warning.) before the ]
% unless the warning was issued with verbose off; long lines are wrapped at
% the window width (also inside that line), and the id may be a hyperlink.
% A warning without that line has an empty identifier. The last warning takes
% lastwarn's identifier and exact message when its identifier is lastwarn's or
% its message is lastwarn's but for line breaks.


if nargin < 3
    lastMsg = '';
    lastId = '';
end
bs = char(8);
w = {};
blocks = regexp(txt, ['\[', bs, 'Warning: (.*?)\]', bs], 'tokens');
if isempty(blocks) && contains(txt, '[Warning: ')
    blocks = regexp(txt, '\[Warning: (.*?)\]\s*(?:\n|$)', 'tokens');
end
typeLine = ['^(.*?)\s*\(Type\s+"(?:<a href="[^"]*">)?warning\s+off\s+([^"<]*?)(?:</a>)?"\s+', ...
    'to\s+suppress\s+this\s+warning\.\)\s*$'];
for i = 1 : numel(blocks)
    body = blocks{i}{1};
    t = regexp(body, typeLine, 'tokens', 'once');
    if isempty(t)
        w{end + 1} = struct('identifier', '', 'message', strtrim(body)); %#ok<AGROW>
    else
        w{end + 1} = struct('identifier', strtrim(t{2}), 'message', t{1}); %#ok<AGROW>
    end
end
if ~isempty(w)
    ident = w{end}.identifier;
    same = strcmp(regexprep(strtrim(w{end}.message), '\s+', ' '), regexprep(strtrim(lastMsg), '\s+', ' '));
    if (strcmp(ident, lastId) && (~isempty(lastId) || same)) || (isempty(ident) && same)
        if isempty(lastId)
            lastId = ident;
        end
        w{end}.identifier = lastId;
        w{end}.message = lastMsg;
    end
end
n = numel(strfind(txt, '[Warning: ')) + numel(strfind(txt, ['[', bs, 'Warning: ']));
if n > numel(blocks)
    w{end + 1} = struct('identifier', 'pkpy:oracle:unparsedWarning', ...
        'message', sprintf('%d warnings could not be parsed; see stdout.txt', n - numel(blocks)));
end

end
