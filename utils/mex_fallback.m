function [records] = mex_fallback(mexName, ME, designedIDs)
% record that the MATLAB code runs instead of the mex function mexName,
% after the mex threw ME, or without calling it (ME empty, e.g., where the
% code checks that the mex exists). A mex that exists for this platform warns
% (identifier PetaKit5D:mexFallback) once per mex, error identifier and
% MATLAB session, since its failure (wrong arguments, a library that doesn't
% load, ...) otherwise only makes the processing slower (per identifier, so
% that a failure that doesn't matter, e.g., a missing file, doesn't hide a
% later one that does). No warning for a mex that isn't built for
% this platform, or a rejection by design: an error identifier in designedIDs
% (e.g., a data type the mex doesn't support).
%
% Every call is recorded by mex, kind ('failed': the mex exists, 'missing':
% it doesn't, 'designed') and error identifier, with the first error message
% and a count. Tests read and reset the records:
%   records = mex_fallback('-log')   struct array (mex, kind, identifier, message, count)
%   mex_fallback('-reset')           clear the records and the warnings given
%
% mexName: char, the mex function, e.g. 'crop3d_mex'.
% ME: MException or [] (default: []).
% designedIDs: cell array of char (default: {}). Error identifiers of rejections by design.


persistent fallbackLog warned

if ~isstruct(fallbackLog)
    fallbackLog = struct('mex', {}, 'kind', {}, 'identifier', {}, 'message', {}, 'count', {});
    warned = {};
end

switch mexName
    case '-log'
        records = fallbackLog;
        return;
    case '-reset'
        fallbackLog = struct('mex', {}, 'kind', {}, 'identifier', {}, 'message', {}, 'count', {});
        warned = {};
        return;
end

if nargin < 2
    ME = [];
end
if nargin < 3
    designedIDs = {};
end

identifier = '';
message = '';
if ~isempty(ME)
    identifier = ME.identifier;
    message = strtrim(regexprep(ME.message, '\s*\n\s*', ' '));
end

if ~isempty(identifier) && any(strcmp(identifier, designedIDs))
    kind = 'designed';
elseif mex_exists(mexName)
    kind = 'failed';
else
    kind = 'missing';
end

ind = find(strcmp({fallbackLog.mex}, mexName) & strcmp({fallbackLog.kind}, kind) & strcmp({fallbackLog.identifier}, identifier), 1);
if isempty(ind)
    fallbackLog(end + 1) = struct('mex', mexName, 'kind', kind, 'identifier', identifier, 'message', message, 'count', 1);
else
    fallbackLog(ind).count = fallbackLog(ind).count + 1;
end

% a warning counts as given only if it shows (not where the caller turned
% warnings off)
key = [mexName, '|', identifier];
if ~strcmp(kind, 'failed') || any(strcmp(warned, key))
    return;
end
ws = warning('query', 'PetaKit5D:mexFallback');
if strcmp(ws.state, 'off')
    return;
end
warned{end + 1} = key;

% the function that called the mex, for the warning
caller = '';
try
    st = dbstack(1);
    if ~isempty(st)
        caller = sprintf(' in %s (line %d)', st(1).name, st(1).line);
    end
catch
end
if isempty(ME)
    reason = sprintf('was not used%s', caller);
else
    reason = sprintf('failed%s: %s', caller, message);
end
% one line, without the backtrace
bt = warning('off', 'backtrace');
warning('PetaKit5D:mexFallback', '%s %s; the MATLAB code runs instead.', mexName, reason);
warning(bt);

end


function [tf] = mex_exists(mexName)
% a mex file on the path, or a function in a package (a name with dots, e.g.,
% MATLAB's built-in matlab.images.internal.resize.imresizemex)

tf = exist(mexName, 'file') == 3;
if ~tf && any(mexName == '.')
    tf = ~isempty(which(mexName));
end

end
