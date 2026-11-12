function [outs, err, txt, tCall] = call(c, args, work)
% call a case's function (or MATLAB code) and capture its command window
% text, its outputs and its error


nout = c.nargout;
outs = {};
err = [];
tCall = NaN;
if ~isempty(c.code)
    code = c.code;
    txt = evalc('[outs, err, tCall] = pkpy_oracle.call_code(code, args, work, nout);');
else
    f = c.func;
    txt = evalc('[outs, err, tCall] = pkpy_oracle.call_function(f, args, nout);');
end

end
