function [out] = pkpy_x27_probe(fname, base, names, probes)
% How a function's argument checks take a battery of calls (the port's checks of third-party
% code are made from these results).
%
% For each probe i, calls fname(base{:}, probes{i}{:}) in a folder of its own, with its
% command-window output captured and dropped, and sets out.(names{i}) to 'accepted' when the
% call got past the function's inputParser (whatever followed: a result or an error of the
% function's own), else to the identifier of inputParser's error. Only identifiers are kept,
% never error messages: inputParser's messages print the function's validators.

out = struct();
here = pwd;
for i = 1 : numel(probes)
    d = fullfile(here, sprintf('probe%03d', i));
    mkdir(d);
    cd(d);
    out.(names{i}) = one_call(fname, [base, probes{i}]);
    cd(here);
    close all force
end

end


function [r] = one_call(fname, args)
try
    evalc('feval(fname, args{:});');
    r = 'accepted';
catch err
    if startsWith(err.identifier, 'MATLAB:InputParser:')
        r = err.identifier;
    else
        r = 'accepted';
    end
end
end
