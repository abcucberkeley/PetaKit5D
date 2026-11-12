function [problem] = waitsOfASecondProblem(waits, n)
% what is wrong with the waits of a wrapper between its looks for results or
% inputs, '' if nothing: they must be n waits, the first of a second and none
% of 2 s or more. On a fast machine each is a second; wait_for_files waits ten
% times as long as a look for the files that took more than a tenth of a
% second (e.g., 1.19 s after a look of 0.119 s on a slow CI runner), and the
% last wait of a round is what is left of it. Waits of a whole round between
% looks (3 s or more), which the wrappers used to make, are wrong.
%
% waits: the durations of the pause calls in seconds; n: the number of waits.


problem = '';
desc = mat2str(reshape(waits, 1, []), 4);
if numel(waits) ~= n
    problem = sprintf('There were %d waits instead of %d: %s s.', numel(waits), n, desc);
elseif n > 0 && waits(1) ~= 1
    problem = sprintf('The first wait was not a second: %s s.', desc);
elseif any(waits <= 0 | waits >= 2)
    problem = sprintf('A wait was not between 0 and 2 s: %s s.', desc);
end

end
