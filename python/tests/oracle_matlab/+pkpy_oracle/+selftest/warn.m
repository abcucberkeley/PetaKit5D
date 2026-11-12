function [y] = warn(x)
% issue two warnings with identifiers and return 2 x


warning('pkpy:selftest:first', 'first warning');
warning('pkpy:selftest:second', 'second warning: %d values', numel(x));
y = x * 2;

end
