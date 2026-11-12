function [x] = fail(n)
% warn, then raise an error with an identifier


warning('pkpy:selftest:beforeFail', 'about to fail with %d', n);
error('pkpy:selftest:fail', 'failed on purpose with %d', n);
x = n; %#ok<UNRCH>

end
