function [out] = pkpy_parse_probe(fname, base, calib, names, probes)
% How a function's inputParser takes a battery of calls, with the errors that its checks raise
% themselves told apart (the port's checks of third-party code are made from these results).
%
% For each probe i, calls fname(base{:}, probes{i}{:}) in a folder of its own, with its
% command-window output captured and dropped, and sets out.(names{i}) to:
%   'accepted': no error was raised while the function's inputParser parsed the call (the
%       call returned, or failed before or after the parse);
%   the identifier of an error of inputParser's own ('MATLAB:InputParser:...');
%   'raised(<identifier>)': an error that a check raised itself, which inputParser passes on
%       with the check's identifier ('raised()' for one without an identifier).
%
% inputParser raises every error of the parse, its own and its checks', as an error of its
% caller at the call of parse: the error's stack starts at that call, with no frame of
% inputParser or of the check. The call fname(base{:}, calib{:}) must fail there with an error
% of inputParser's own (a parameter's name without its value, say): the first frame of its
% stack is the parse site. An error was raised in the parse when its stack starts at that
% frame (file, function and line) and has as many frames. Only identifiers and the result of
% that comparison are kept, never a message (inputParser's messages print the validator) or a
% frame's name (an anonymous validator's text).

out = struct();
here = pwd;
ref = one_call(fname, [base, calib], fullfile(here, 'probe000'), here);
if isempty(ref) || ~own_error(ref)
    error('pkpy:probe:noParseSite', ...
        'The calibration call of %s raised no error of inputParser''s own.', fname);
end
for i = 1 : numel(probes)
    err = one_call(fname, [base, probes{i}], fullfile(here, sprintf('probe%03d', i)), here);
    if isempty(err) || ~same_site(err, ref)
        out.(names{i}) = 'accepted';
    elseif own_error(err)
        out.(names{i}) = err.identifier;
    else
        out.(names{i}) = ['raised(' err.identifier ')'];
    end
end

end


function [err] = one_call(fname, args, d, here)
% the call's error ([] if none), made in folder d with its output captured and dropped
mkdir(d);
cd(d);
err = [];
try
    evalc('feval(fname, args{:});');
catch e
    err = e;
end
cd(here);
close all force
end


function [tf] = own_error(err)
tf = startsWith(err.identifier, 'MATLAB:InputParser:');
end


function [tf] = same_site(a, b)
% whether two errors were raised at the same place: the first frames of their stacks and the
% number of frames
tf = ~isempty(a.stack) && numel(a.stack) == numel(b.stack) ...
    && strcmp(a.stack(1).file, b.stack(1).file) && strcmp(a.stack(1).name, b.stack(1).name) ...
    && a.stack(1).line == b.stack(1).line;
end
