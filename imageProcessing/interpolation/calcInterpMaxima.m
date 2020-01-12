%[fmax, xmax] = calcInterpMaxima(s, varargin) returns the local maxima of the cubic spline interpolation of s
%
% Examples:
% 
% s = rand(1,10);
% calcInterpMaxima(s, 0, 'Display', true);
% 
% s = rand(4,8);
% [fmax, xmax] = calcInterpMaxima(s, 0.0, 'Display', true);

% Francois Aguet, 07/21/2013

function [fmax, xmax, c] = calcInterpMaxima(f0, varargin)

ip = inputParser;
ip.addOptional('lambda', 0, @(x) isscalar(x) && x>=0);
ip.addParamValue('Display', false, @islogical);
ip.parse(varargin{:});

dims = size(f0);
dims(dims==1) = [];
switch numel(dims)
    case 1
        nx = numel(f0);
        
        % compute coefficients
        c = computeBSplineCoefficients(f0, ip.Results.lambda);
        cx = [c(2) c c(end-1)];
        
        % for each interval, calculate solution
        xmax = cell(1,nx-1);
        for i = 1:nx-1
            xmax{i} = i+getSol(cx(i+(0:3)));
        end
        xmax = [xmax{:}];
        fmax = arrayfun(@(i) interpBSplineValue(i, c, 3, 'symmetric'), xmax);
        
    case 2
        % compute coefficients
        c = computeBSplineCoefficients(f0, ip.Results.lambda);
        [ny,nx] = size(f0);
        
        % pad coefficients for mirroring
        cx = padarrayXT(c, [1 1], 'symmetric');
        
        %opts2 = optimoptions('fmincon', 'Display', 'off',...
        %    'MaxFunEvals', 1e4, 'MaxIter', 1e4, ...
        %    'TolX', 1e-6, 'TolFun', 1e-6, 'TolCon', 1e-6,...
        %    'Algorithm', 'interior-point', 'GradObj', 'on', 'GradConstr', 'on');
                
        % loop over grid, and calculate maximum for each set of points
        xmax = cell(ny-1,nx-1);
        ymax = cell(ny-1,nx-1);
        fmax = cell(ny-1,nx-1);
        for i = 1:nx-1
            for j = 1:ny-1
                %
                %i = 2;
                %j = 1;
                
                C0 = cx(j+(0:3),i+(0:3));
                [dx,dy,val] = iterateSol(C0);
                if ~isempty(dx)                   
                    xmax{j,i} = i+dx;
                    ymax{j,i} = j+dy;
                    fmax{j,i} = val;
                end
                
                %[dx,f,exitflag] = fmincon(@(dx) cost2D(dx, C0), 0.51,...
                %    [], [], [], [], 0, 1, @(dx) constr2D(dx, C0), opts2);
                %if (exitflag==2 || exitflag==1) % valid solution found
                %    cy = C0*b3(dx)';
                %    dy = getSol(cy(:)');
                %    if ~isempty(dy) %&& dy~=0 &&dy~=1
                %        % calculate Hessian
                %        fxx = b3(dy)*C0*d2b3(dx)';
                %        fyy = d2b3(dy)*C0*b3(dx)';
                %        fxy = db3(dy)*C0*db3(dx)';
                %        D = fxx*fyy-fxy^2;
                %        if D>0 && fxx<0 % local max.
                %            xmax{j,i} = i+dx;
                %            ymax{j,i} = j+dy;
                %            fmax{j,i} = -f;
                %        end
                %    end
                %end
            end
        end
        xmax = [vertcat(xmax{:}) vertcat(ymax{:})];
        fmax = vertcat(fmax{:});
        % discard maxima at borders
        if ~isempty(xmax)
            d = 1e-10;
            bidx = xmax(:,1)<=1+d | xmax(:,1)>=nx-d | xmax(:,2)<=1+d | xmax(:,2)>=ny-d;
            xmax(bidx,:) = [];
            fmax(bidx,:) = [];
        end
end

if ip.Results.Display
    switch numel(dims)
        case 1
            
            figure;
            hold on;
            plot([0 nx+1], [0 0], 'k--', 'HandleVisibility', 'off');
            plot(f0, 'k.');
            
            % value of the 1st derivative on grid
            v0 = conv([c(2) c c(end-1)], [0.5 0 -0.5], 'valid');
            plot(v0, 'b.', 'HandleVisibility', 'off');
            
            xi = 1:0.01:nx;
            [fi, dfi] = binterp(f0, xi);
            plot(xi, fi, 'k-');
            plot(xi, dfi, 'b-');
            plot(xmax, fmax, 'ro', 'MarkerSize', 10);
            YLim = get(gca, 'YLim');
            set(gca, 'XLim', [0.5 nx+0.5], 'YLim', max(abs(YLim))*[-1 1]);
            legend(' Samples', ' Cubic spline interpolation', ' Derivative', 'Location', 'NorthOutside');
            
        case 2
            % upsample 100x
            xa = 1:0.01:nx;
            ya = 1:0.01:ny;
            [xi,yi] = meshgrid(xa,ya);
            if ip.Results.lambda==0
                Si = binterp(f0, xi, yi); % use C version for speed
            else
                Si = xinterp(c, 10);
            end
            figure; imagesc(xa, ya, Si); colormap(jet(256)); axis image; colorbar;
            hold on;
            if ~isempty(xmax)
                plot(xmax(:,1), xmax(:,2), 'kx', 'MarkerSize', 10);
            end
    end
end



function sol = getSol(c)
A = c(1)-3*c(2)+3*c(3)-c(4);
B = c(1)-2*c(2)+c(3);
delta = -c(1)*c(2) + 4*c(2)^2 - 7*c(2)*c(3) + 4*c(3)^2 + c(1)*c(4) - c(3)*c(4);
if delta<0
    sol = [];
else
    delta = sqrt(delta);
    sol = [B-delta B+delta]/A;
    sol(abs(sol)<eps) = 0;
    
    % check 2nd derivative
    d2 = arrayfun(@(i) sum(c.*d2b3(i)), sol);
    sol = sol(d2<0);
    
    % check bounds
    sol = sol(sol>=0 & sol<=1);
end


function [x, y, val] = iterateSol(C)

% Convergence criteria: change in dx and dy less than Tol
tol = 1e-10;
maxIter = 20;

% Start with initial values at both interval boundaries since
% two solutions for dF/dx==0 may exist within interval
yv = [0.1 0.99 0.5];
ns = numel(yv);

D = NaN(ns,2);
for d = 1:ns
    maxDelta = 1;
    iter = 1;

    % assuming dy, calculate dx first
    y0 = yv(d);
    cx = b3(y0)*C;
    x0 = getSol(cx);
    if ~isempty(x0) % iterate
        x = x0;
        while maxDelta>tol && iter<=maxIter
        
            cy = C*b3(x)';
            y = getSol(cy(:)');
            if isempty(y)
                break;
            end
            dy = y-y0;
            y0 = y;
        
            if length(y)>1 && y(1)==y(2)
                y = y(1);
            end
            cx = b3(y)*C;
            x = getSol(cx);
            if isempty(x)
                break;
            end
            dx = x-x0;
            x0 = x;
            
            maxDelta = max(abs([dx dy]));
            iter = iter + 1;
            
        end
        if ~isempty(x) && ~isempty(y)
            D(d,:) = [x y];
        end        
    end
end
D(isnan(D(:,1)),:) = [];
% evaluate valid solutions

ns = size(D,1);
val = NaN(ns,1);
for i = 1:ns
    % calculate Hessian
    fxx = b3(D(i,2))*C*d2b3(D(i,1))';
    fyy = d2b3(D(i,2))*C*b3(D(i,1))';
    fxy = db3(D(i,2))*C*db3(D(i,1))';
    d = fxx*fyy-fxy^2;
    if d>0 && fxx<0 % local max.
        val(i) = b3(D(i,2))*C*b3(D(i,1))';
    end
end
% idx = ~isnan(val);
% val = val(idx);
[val,idx] = max(val);
x = D(idx,1);
y = D(idx,2);



% cubic spline kernel (reversed for convolution!)
function k = b3(dx)
k = [(1-dx)^3 4+3*(dx-2)*dx^2 1+3*dx*(1+dx-dx^2) dx^3]/6;

function k = db3(dx)
k = [-0.5*(dx-1)^2 0.5*dx*(3*dx-4) 0.5+dx-1.5*dx^2 0.5*dx^2];

function k = d2b3(dx)
k = [1-dx -2+3*dx 1-3*dx dx];



function [v, g] = cost2D(dx, C)
% Given dx, analytically solve for dy
cy = C*b3(dx)';
dy = getSol(cy(:)');

% if no maxima are found, choose highest boundary
if isempty(dy)
    if sum(b3(0).*cy') > sum(b3(1).*cy')
        dy = 0;
    else
        dy = 1;
    end
end

v = -b3(dy)*cy; % negative sign for minization
g = -b3(dy)*C*db3(dx)';

function [c, ceq, GC, GCeq] = constr2D(dx, C)
c = [];
GC = [];

cy = C*b3(dx)';
dy = getSol(cy(:)');

% if no maxima are found, choose highest boundary
if isempty(dy)
    if sum(b3(0).*cy') > sum(b3(1).*cy')
        dy = 0;
    else
        dy = 1;
    end
end

% force gradient to zero
ceq = -b3(dy)*C*db3(dx)';
GCeq = -b3(dy)*C*d2b3(dx)';



function g = xinterp(c, f)
[ny,nx] = size(c);
c = padarrayXT(c, [2 2], 'symmetric');

df = 1/f;
x = 1:df:nx;
y = 1:df:ny;
nx = numel(x);
ny = numel(y);
g = zeros(ny,nx);

for i = 1:nx
    for j = 1:ny
        x0 = floor(x(i));
        dx = x(i)-x0;
        y0 = floor(y(j));
        dy = y(j)-y0;
        c0 = c(y0+(1:4), x0+(1:4));
        g(j,i) = b3(dy)*c0*b3(dx)';
    end
end
