function s = pkpy_f02_fileparts(p)
% fileparts' three outputs as 'path|name|ext'
[a, b, c] = fileparts(p);
s = sprintf('%s|%s|%s', a, b, c);
end
