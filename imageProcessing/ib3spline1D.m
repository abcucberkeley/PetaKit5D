function ima = ib3spline1D(coeffs, nx)

if nargin == 1 || isempty(nx)
    nx = size(coeffs, 2);
end

[nyCoeffs, nxCoeffs] = size(coeffs);

% Rescale coordinates
xScaled = (1:nx) / (nx / nxCoeffs);

% Compute the interpolation indexes
xIndex = bsxfun(@plus, floor(xScaled), (-1:2)');

% Compute weights
t = xScaled - xIndex(2,:);

w = zeros(size(xIndex));
w(4,:) = (1/6) * t.^3;
w(1,:) = (1/6) + (1/2) * t .* (t - 1) - w(4,:);
w(3,:) = t + w(1,:) - 2 * w(4,:);
w(2,:) = 1 - w(1,:) - w(3,:) - w(4,:);

% Add mirror condition at border 
coeffsPadded = padarrayXT(coeffs, [0 2], 'symmetric', 'both');

ima = zeros(nyCoeffs,nx);

for k = 1:4
    ima = ima + repmat(w(k,:), nyCoeffs,1) .* coeffsPadded(:, xIndex(k,:) + 2);
end
