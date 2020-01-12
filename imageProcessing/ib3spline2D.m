function ima = ib3spline2D(coeffs, imaSize)

if nargin == 1 || isempty(imaSize)
    [ny,nx] = size(coeffs);
else
    ny = imaSize(1);
    nx = imaSize(2);
end

ima = ib3spline1D(coeffs,nx);
ima = ib3spline1D(ima',ny)';
