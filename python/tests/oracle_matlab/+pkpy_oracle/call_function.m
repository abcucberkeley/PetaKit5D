function [outs, err, tCall] = call_function(f, args, nout)
% feval with nout outputs; an error is returned, not raised


outs = cell(1, nout);
err = [];
t = tic;
try
    if nout > 0
        [outs{1 : nout}] = feval(f, args{:});
    else
        feval(f, args{:});
    end
catch ME
    err = ME;
end
tCall = toc(t);

end
