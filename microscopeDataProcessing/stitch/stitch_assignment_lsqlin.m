function [x] = stitch_assignment_lsqlin(C, d, A, b, Aeq, beq)
% the shifts of the tiles along one axis for the global shift assignment: the
% minimum of ||C x - d||^2 subject to A x <= b and Aeq x = beq (as lsqlin),
% whose rows are weighted differences of the shifts of two tiles, made unique:
%
% 1. Determined values are the exact optimum: the least squares of C alone if
%    it holds the inequalities (with the groups placed as in 2.), otherwise
%    lsqlin's solution refined on its active set (the active rows of A and the
%    rows of Aeq as equalities).
% 2. Groups of tiles that C and Aeq don't connect (tiles whose pairs are all
%    below the xcorr threshold, or groups linked only by such pairs) are free
%    in the objective. They keep their stage positions relative to each other
%    as far as the bounds allow: their offsets minimize the sum of squares of
%    the rows of A between groups (the shift differences of their overlapping
%    pairs) within the bounds of these rows. So a free tile gets the mean of
%    the shifts of the tiles it overlaps (each pair alike), a free group the
%    mean over its pairs to other groups of the neighbor's shift less its
%    tile's, and adjacent free groups are placed together.
% 3. Each group of overlapping tiles (connected by the rows of C, A and Aeq)
%    has mean zero.
%
% Without rows in C (every pair below the threshold) the tiles keep their stage
% positions (x = 0). It runs on one thread, so that the result does not depend
% on MATLAB's thread count.

if nargin < 5
    Aeq = [];
    beq = [];
end
n = max([size(C, 2), size(A, 2), size(Aeq, 2)]);
if n == 0
    x = zeros(0, 1);
    return;
end
C = sparse(C);
A = sparse(A);
if isempty(Aeq)
    Aeq = sparse(0, n);
    beq = zeros(0, 1);
end
Aeq = sparse(Aeq);
d = double(full(d(:)));
b = double(full(b(:)));
beq = double(full(beq(:)));

nt = maxNumCompThreads(1);
restoreThreads = onCleanup(@() maxNumCompThreads(nt));

% groups of tiles: of the rows of C and Aeq (their offsets are free in the
% objective), and of all rows (a constant added to such a group changes nothing)
ec = row_ends(C);
[ea, ra] = row_ends(A);
ee = row_ends(Aeq);
grp = conncomp(graph([ec(:, 1); ee(:, 1)], [ec(:, 2); ee(:, 2)], [], n))';
lin = conncomp(graph([ec(:, 1); ea(:, 1); ee(:, 1)], [ec(:, 2); ea(:, 2); ee(:, 2)], [], n))';
tol = 1e-9 * (1 + max([0; abs(b)]));

x = zeros(n, 1);
ok = false;
if isempty(ec)
    % no rows in C: the stage positions
    ok = true;
elseif isempty(ee)
    % the optimum of C alone (the first tile of each group fixed)
    keep = true(n, 1);
    keep(first_members(grp)) = false;
    Ck = C(:, keep);
    x(keep) = (Ck' * Ck) \ (Ck' * d);
    [x, ok] = place_groups(x, A, b, ea, ra, grp, lin, tol);
    ok = ok && all(A * x <= b + tol);
end
if ~ok
    % lsqlin, refined on its active set (the first tile of each group of
    % overlapping tiles fixed)
    keep = true(n, 1);
    keep(first_members(lin)) = false;
    [u, ~, ~, ~, ~, lambda] = lsqlin(C(:, keep), d, A(:, keep), b, Aeq(:, keep), beq);
    u = refine_active_set(C(:, keep), d, A(:, keep), b, Aeq(:, keep), beq, u, lambda);
    x = zeros(n, 1);
    x(keep) = u;
    x = place_groups(x, A, b, ea, ra, grp, lin, tol);
end

% mean zero in each group of overlapping tiles
gm = accumarray(lin, x) ./ accumarray(lin, 1);
x = x - gm(lin);

end


function [x, ok] = place_groups(x, A, b, ea, ra, grp, lin, tol)
% the offsets of the groups of C and Aeq (the group of the first tile of each
% group of overlapping tiles fixed) with the least sum of squares of the rows of
% A between groups, within the bounds of these rows; ok is false if no offsets
% hold these bounds
ok = true;
cross = grp(ea(:, 1)) ~= grp(ea(:, 2));
k = max(grp);
mov = true(k, 1);
mov(grp(first_members(lin))) = false;
if ~any(cross) || ~any(mov)
    return;
end
rows = ra(cross);
P = sparse(1 : numel(grp), grp, 1, numel(grp), k);
G = A(rows, :) * P(:, mov);
r0 = A(rows, :) * x;
s0 = b(rows) - r0;
c = -((G' * G) \ (G' * r0));
if any(G * c > s0 + tol)
    % a bound holds them: lsqlin, refined on its active set
    opts = optimoptions('lsqlin', 'Display', 'off');
    [c, ~, ~, flag, ~, lambda] = lsqlin(G, -r0, G, s0, [], [], [], [], [], opts);
    if flag <= 0
        ok = false;
        return;
    end
    c = refine_active_set(G, -r0, G, s0, sparse(0, size(G, 2)), zeros(0, 1), c, lambda);
end
cm = zeros(k, 1);
cm(mov) = c;
x = x + cm(grp);

end


function [u] = refine_active_set(C, d, G, h, E, e, u, lambda)
% lsqlin's solution u refined on its active set (the rows whose multiplier
% exceeds their slack); as it is without multipliers (an infeasible problem)
if ~isstruct(lambda)
    return;
end
z = lambda.ineqlin;
s = h - G * u;
act = z > max(s, 1e-9 * (1 + max(abs(z), [], 'all')));
u = polish(C, d, G, h, E, e, u, act, z, lambda.eqlin);

end


function [u, act] = polish(C, d, G, h, E, e, u, act, z, y)
% the optimum on lsqlin's active set: the least squares of C with the active
% rows of G and the rows of E as equalities, by iterative refinement of the
% regularized KKT system; rows with a negative multiplier leave, violated rows
% join; lsqlin's solution if that does not end optimal
r = size(C, 2);
H = C' * C;
f = -(C' * d);
rho = 1e-9 * max(1, full(max(abs(H), [], 'all')));
tol = 1e-9 * (1 + max([0; abs(h)]));
zw = z;
act0 = act;
for pass = 1 : 30
    B = [G(act, :); E];
    c = [h(act); e];
    w = [max(zw(act), 0); y];
    k = size(B, 1);
    % (quasi-definite; ill-conditioned by design, which the refinement handles)
    K = decomposition([H + rho * speye(r), B'; B, -rho * speye(k)], 'lu', 'CheckCondition', false);
    v = u;
    for it = 1 : 6
        sol = K \ [-(H * v + f + B' * w); c - B * v];
        v = v + sol(1 : r);
        w = w + sol(r + 1 : end);
    end
    if ~all(isfinite(v))
        break;
    end
    na = nnz(act);
    lam = w(1 : na);
    slack = h - G * v;
    neg = lam < -tol;
    viol = slack < -tol & ~act;
    if ~any(neg) && ~any(viol)
        u = v;
        return;
    end
    idx = find(act);
    act(idx(neg)) = false;
    zw(idx) = lam .* ~neg;
    if any(viol)
        % the most violated row joins (with a multiplier of 0)
        slack(~viol) = inf;
        [~, j] = min(slack);
        act(j) = true;
        zw(j) = 0;
    end
end
act = act0;

end


function [e, r] = row_ends(M)
% the two tiles (columns) of each row with entries, and the rows' indices
[ri, ci] = find(M);
ri = ri(:);
ci = ci(:);
[~, o] = sortrows([ri, ci]);
ri = ri(o);
ci = ci(o);
cnt = accumarray(ri, 1, [size(M, 1), 1]);
if any(cnt ~= 0 & cnt ~= 2)
    error('stitch_assignment_lsqlin:rowError', 'Each row must be a difference of the shifts of two tiles.');
end
e = reshape(ci, 2, [])';
r = ri(1 : 2 : end);
if isempty(e)
    e = zeros(0, 2);
    r = zeros(0, 1);
end

end


function [f] = first_members(g)
% the first (lowest) index of each group
[~, f] = unique(g, 'first');

end
