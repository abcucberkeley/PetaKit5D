%[c] = computeBSplineCoefficients(s, varargin) computes spline coefficients for input s.
% 
% Inputs: 
%           s: input signal, 1D or 2D
%
% Options:
%      lambda: regularization parameter for smoothing splines.
%              Smoothing-spline coefficients are calculated in the Fourier domain.
%    'Degree': 1-3. Selects linear, quadratic, resp. cubic spline interpolation.
%              Smoothing splines are implemented for degree 3 only.
%      'Mode': {'Fourier'}|'Spatial' selects method for coefficient calculation.
%  'Boundary': 'periodic'|{'symmetric'} selects the boundary conditions.
%
% Outputs:
%           c: spline coefficients
%
% Notes: This function is based on the formalism/algorithm described in:
% [1] Unser, IEEE Signal Proc. Mag. 16(6), pp. 22-38, 1999
% [2] Unser et al., IEEE Trans. Signal Proc. 41(2), pp. 834-848, 1993
%     For details, see Section IV. B, Eq. 4.3

% Francois Aguet, 2010 (Last modified: 07/21/2013)

function c = computeBSplineCoefficients(s, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('s', @isnumeric);
ip.addOptional('lambda', 0, @(x) isscalar(x) && x>=0);
ip.addParamValue('Degree', 3, @(x) ismember(x, 1:3));
ip.addParamValue('Mode', 'fourier', @(x) any(strcmpi(x, {'fourier', 'spatial'})));
ip.addParamValue('Boundary', 'symmetric', @(x) any(strcmpi(x, {'symmetric', 'periodic', 'replicate', 'zeros'})));
ip.parse(s, varargin{:});

lambda = ip.Results.lambda;
mode = ip.Results.Mode;
if lambda~=0
    mode = 'fourier';
end

if strcmpi(mode, 'fourier')
    dims = size(s);
    dims(dims==1) = [];
    switch numel(dims)
        case 1
            N = numel(s);
            
            % mirror signal
            if strcmpi(ip.Results.Boundary, 'symmetric')
                s = [s s(N-1:-1:2)];
                M = 2*N-2;
            else % periodic
                M = N;
            end
            
            % frequency vector
            w = (0:M-1)*2*pi/M;
            
            % Fourier transform of spatial input signal
            S = fft(s);
            
            % Smoothing spline pre-filter
            H = 3 ./ (2+cos(w)+6*lambda*(cos(2*w)-4*cos(w)+3));
            
            % Spline coefficients
            c = real(ifft(S.*H));
            c = c(1:N);
        case 2
            [ny,nx] = size(s);
            
            % mirror signal
            if strcmpi(ip.Results.Boundary, 'symmetric')
                s = [s s(:,nx-1:-1:2)];
                s = [s; s(ny-1:-1:2,:)];
                mx = 2*nx-2;
                my = 2*ny-2;
            else % periodic
                mx = nx;
                my = ny;
            end
            
            % frequency vectors
            wx = (0:mx-1)*2*pi/mx;
            wy = (0:my-1)'*2*pi/my;
            
            % Fourier transform of spatial input signal
            S = fft2(s);
            
            % Smoothing spline pre-filters
            Hx = 3 ./ (2+cos(wx)+6*lambda*(cos(2*wx)-4*cos(wx)+3));
            Hy = 3 ./ (2+cos(wy)+6*lambda*(cos(2*wy)-4*cos(wy)+3));
            
            % Spline coefficients
            c = real(ifft2((Hy*Hx).*S));
            c = c(1:ny,1:nx);
    end
else % calculate in spatial domain
    if ismember(ip.Results.Degree, [2 3])
        [ny,nx] = size(s);
        switch ip.Results.Degree
            case 3
                z1 = -2+sqrt(3);
                c0 = 6;
            case 2
                z1 = -3+2*sqrt(2);
                c0 = 8;
        end
        
        % Recursively compute coefficients
        if nx>1
            cp = zeros(ny,nx);
            cn = zeros(ny,nx);
            
            cp(:,1) = getCausalInitValue(s, z1, ip.Results.Boundary);
            for k = 2:nx
                cp(:,k) = s(:,k) + z1*cp(:,k-1);
            end
            cn(:,nx) = getAntiCausalInitValue(cp, z1, ip.Results.Boundary);
            for k = nx-1:-1:1
                cn(:,k) = z1*(cn(:,k+1)-cp(:,k));
            end
            c = c0*cn;
        else
            c = s;
        end
        if ny>1
            c = c'; % quick hack: transpose
            cp = zeros(nx,ny);
            cn = zeros(nx,ny);
            
            cp(:,1) = getCausalInitValue(c, z1, ip.Results.Boundary);
            for k = 2:ny
                cp(:,k) = c(:,k) + z1*cp(:,k-1);
            end
            cn(:,ny) = getAntiCausalInitValue(cp, z1, ip.Results.Boundary);
            for k = ny-1:-1:1
                cn(:,k) = z1*(cn(:,k+1)-cp(:,k));
            end
            c = c0*cn'; % transpose back
        end
    else % n = 1
        c = s;
    end
end



function c0 = getCausalInitValue(s, a, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('s', @isnumeric);
ip.addRequired('a', @isscalar);
ip.addOptional('Boundary', 'symmetric', @(x) any(strcmpi(x, {'symmetric', 'periodic', 'replicate', 'zeros'})));
ip.parse(s, a, varargin{:});
N = size(s,2);

switch ip.Results.Boundary
    case 'symmetric'
        k = repmat(0:2*N-3, [size(s,1) 1]);
        s = [s s(:,end-1:-1:2)];
        c0 = sum(s.*a.^k, 2) / (1 - a^(2*N-2));
    case 'periodic'
        k = repmat(0:N-1, [size(s,1) 1]);
        s = [s(:,1) s(:,N:-1:2)];
        c0 = sum(s.*a.^k, 2) / (1-a^N);
    case 'zeros'
        c0 = s(:,1);
    case 'replicate'
        c0 = s(:,1)/(1-a);
end



function c0 = getAntiCausalInitValue(c, a, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('c', @isnumeric);
ip.addRequired('a', @isscalar);
ip.addOptional('Boundary', 'symmetric', @(x) any(strcmpi(x, {'symmetric', 'periodic', 'replicate', 'zeros'})));
ip.parse(c, a, varargin{:});

N = size(c,2);

switch ip.Results.Boundary
    case 'symmetric'
        c0 = (a/(a*a-1))*(c(:,N)+a*c(:,N-1)); % [1] Box 2 (has errors)
        % Equivalent: (needs input signal), see [2], eqs. 2.3 & 2.5
        %c0 = a/(a*a-1) * (2*c(:,N)-s(N));
    case 'periodic'
        k = repmat(0:N-1, [size(c,1) 1]);
        c0 = -a/(1-a^N) * sum(a.^k .* c(:,[N 1:N-1]),2);
    case 'zeros'
        c0 = -a*c(:,end);
    case 'replicate'
        c0 = -c(:,end)*a/(1-a);
end
