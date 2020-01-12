% [c] = b3spline1D(img, boundary)
% 
% Computes the 1D cubic B-spline coefficients.
%
% Inputs: 
%           img      : input image
%           boundary : boundary conditions: 'mirror' (default) or 'periodic'

% Francois Aguet, June 2010

function cn = b3spline1D(img, boundary)

if nargin<2
    boundary = 'mirror';
end

% cubic spline parameters
c0 = 6;
z1 = -2+sqrt(3);

N = size(img,2);
cp = zeros(size(img));
cn = zeros(size(img));

if (strcmp(boundary, 'mirror')==1)
    cp(:,1) = getCausalInit_Mirror(img, z1);
    for k = 2:N
        cp(:,k) = img(:,k) + z1*cp(:,k-1);
    end;
    cn(:,N) = getAntiCausalInit_Mirror(cp, z1);
    for k = N-1:-1:1
        cn(:,k) = z1*(cn(:,k+1) - cp(:,k));
    end;
elseif (strcmp(boundary, 'periodic')==1)
    cp(:,1) = getCausalInit_Periodic(img, z1);
    for k = 2:N
        cp(:,k) = img(:,k) + z1*cp(:,k-1);
    end;
    cn(:,N) = getAntiCausalInit_Periodic(cp, z1);
    for k = N-1:-1:1
        cn(:,k) = z1*(cn(:,k+1) - cp(:,k));
    end
else
    error('Boundary: unknown boundary conditions');
end;
cn = c0*cn;


function c0 = getAntiCausalInit_Mirror(img, a)
N = size(img,2);
c0 = (a/(a*a-1))*(img(:,N)+a*img(:,N-1));


function c0 = getAntiCausalInit_Periodic(img, a)
N = size(img,2);
img = [img(:,N) img img(:,1:N-3)];
k = repmat(0:2*N-3, [size(img,1) 1]);
c0 = -a/(1-a^N)*sum(a.^k.*img, 2);


function out = getCausalInit_Mirror(img, a)
N = size(img,2);
k = repmat(0:2*N-3, [size(img,1) 1]);
img = [img img(:,end-1:-1:2)];
out = sum(img.*a.^k,2) / (1 - a^(2*N-2));


function out = getCausalInit_Periodic(img, a)
N = size(img,2);
k = repmat(0:N-1, [size(img,1) 1]);
img = [img(:,1) img(:,N:-1:2)];
out = sum(img.*a.^k,2)/(1-a^N);