function [p] = fullPath(p)
% the full path of file or folder p (relative to the current folder), with
% . and .. and symbolic links resolved and no trailing slash

f = java.io.File(p);
if ~f.isAbsolute()
    % (Java resolves relative paths from the folder MATLAB started in)
    f = java.io.File(pwd, p);
end
p = char(f.getCanonicalPath());

end
