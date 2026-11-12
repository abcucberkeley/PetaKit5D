function [outs, err, tCall] = call_code(pkpyCode__, in, work, pkpyNout__) %#ok<INUSL>
% evaluate a case's MATLAB code; it reads in (cell of the arguments) and
% work (the work folder) and sets out (cell of the outputs). The harness's
% own variables have names the code does not use.


out = {};
pkpyErr__ = [];
pkpyT__ = tic;
try
    eval(pkpyCode__);
catch pkpyME__
    pkpyErr__ = pkpyME__;
end
tCall = toc(pkpyT__);
err = pkpyErr__;
if ~iscell(out)
    out = {out};
end
outs = reshape(out, 1, []);
if isempty(err) && numel(outs) ~= pkpyNout__
    error('pkpy:oracle:nargout', 'the case code set %d outputs, nargout is %d', numel(outs), pkpyNout__);
end

end
