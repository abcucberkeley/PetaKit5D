%v = interpBSplineValue(x, c, n, boundary)
%
% Inputs:
%            x : interpolation coordinates
%            c : function to interpolate
%            n : degree of the spline, 1-3. Default: 3
%     boundary : boundary conditions: 'symmetric' (default) or 'periodic'
%
% Outputs:
%            v : interpolated value at x

% Francois Aguet, 2009 (Last modified: 10/05/2011)

function v = interpBSplineValue(x, c, n, boundary)

xi = floor(x);
dx = x-xi;
nx = length(c);

if (n == 1)
    if (strcmp(boundary, 'periodic')==1)
        x0 = periodize(xi, nx);
        x1 = periodize(xi+1, nx);
    else
        x0 = mirror(xi, nx);
        x1 = mirror(xi+1, nx);
    end
    v = dx*c(x1) + (1.0-dx)*c(x0);
else
    if (n == 2)
        wx = getQuadraticSpline(dx);
    else %if (n == 3)
        wx = getCubicSpline(dx);
    end
    if (strcmp(boundary, 'periodic')==1)
        v = sum(wx.*c([periodize(xi-1, nx) periodize(xi, nx) periodize(xi+1, nx) periodize(xi+2, nx)]));
    else
        switch xi
            case 0
                v = sum(wx.*c([3 2 1 2]));
            case 1
                v = sum(wx.*c([2 1 2 3]));
            case nx-1
                v = sum(wx.*c([nx-2 nx-1 nx nx-1]));
            case nx
                v = sum(wx.*c([nx-1 nx nx-1 nx-2]));
            otherwise
                v = sum(wx.*c([xi-1 xi xi+1 xi+2]));
        end
    end
end



function x = mirror(x, nx)
idx = x<1;
x(idx) = 2-x(idx);
idx = x>nx;
x(idx) = 2*nx-x(idx);

% if (x >= 1 && x <= nx)
%     y = x;
% elseif (x < 1)
%     y = 2-x;
% else
%     y = 2*nx-x;
% end


function y = periodize(x, nx)
y = x;
while (y < 1)
    y = y+nx;
end
while (y > nx)
    y = y-nx;
end


% @param	t argument between 0 and 1.
% @return	4 sampled values of the quadratic B-spline
%			(B3[t+1], B3[t], B3[t-1], B3[t-2])
function v = getQuadraticSpline(t)
if (t < 0.0 || t > 1.0)
    error('Argument t for cubic B-spline outside of expected range [0, 1]: %f',t);
end
if (t <= 0.5)
    v(1) = (t-0.5)*(t-0.5)/2.0;
    v(2) = 0.75 - t*t;
    v(3) = 1.0-v(2)-v(1); %(t+0.5)*(t+0.5)/2.0;
    v(4) = 0.0;
else
    v(1) = 0.0;
    v(2) = (t-1.5)*(t-1.5)/2.0;
    v(4) = (t-0.5)*(t-0.5)/2.0;
    v(3) = 1.0-v(4)-v(2);
end


% @param	t argument between 0 and 1.
% @return	4 sampled values of the cubic B-spline
%    		(B3[t+1], B3[t], B3[t-1], B3[t-2]).
function v = getCubicSpline(t)
if (t < 0.0 || t > 1.0)
    error('Argument t for cubic B-spline outside of expected range [0, 1]: %f',t);
end
t1 = 1.0 - t;
t2 = t * t;
v(1) = (t1 * t1 * t1) / 6.0;
v(2) = (2.0 / 3.0) + 0.5 * t2 * (t-2);
v(4) = (t2 * t) / 6.0;
v(3) = 1.0 - v(4) - v(2) - v(1);
