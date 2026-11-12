function [str, strs, map, fh, ip, sp, sa, i64, specials, cs, chm, lg, es, ec] = make_values()
% MATLAB values Python has no literal for, to test how results are read


str = "hello";
strs = ["a", missing, "ccc"];
map = containers.Map({'b', 'a'}, {1, [2, 3]});
fh = @(x) x + 1;
ip = inputParser;
ip.addParameter('alpha', 3);
ip.addParameter('name', 'x');
ip.parse('alpha', 5);
sp = sparse([1, 3], [2, 2], [1.5, -2], 3, 4);
sa = struct('a', {1, 'two'}, 'b', {[], {3}});
i64 = [intmin('int64'), intmax('int64'), int64(-1)];
specials = [0 / 0, -0, Inf, -Inf, 1e-310, NaN];
cs = complex(single(1.5), single(-0));
chm = ['ab'; 'cd'];
lg = logical([1, 0; 0, 1]);
es = struct([]);
ec = cell(0, 3);

end
