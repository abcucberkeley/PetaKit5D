function [w] = parse_warnings(txt)
% the warnings in captured command window text (verbose on, backtrace off):
% a cell array of structs with identifier and message. MATLAB writes a
% warning as [<backspace>Warning: message, a line (Type "warning off <id>" to
% suppress this warning.), and ]<backspace>; the id may be a hyperlink.


bs = char(8);
w = {};
pat = ['\[', bs, '?Warning: (.*?)\n\(Type "(?:<a href="[^"]*">)?warning off ([^"<]*)(?:</a>)?" ', ...
    'to suppress this warning\.\)\s*\]', bs, '?'];
tok = regexp(txt, pat, 'tokens');
for i = 1 : numel(tok)
    w{end + 1} = struct('identifier', strtrim(tok{i}{2}), 'message', tok{i}{1}); %#ok<AGROW>
end
n = numel(regexp(txt, ['\[', bs, '?Warning: '], 'start'));
if n > numel(tok)
    w{end + 1} = struct('identifier', 'pkpy:oracle:unparsedWarning', ...
        'message', sprintf('%d warnings could not be parsed; see stdout.txt', n - numel(tok)));
end

end
