%[fi, fi_dx, fi_d2x] = binterp(f, xi, borderCondition) returns the cubic spline interpolation of 1D or 2D input signal.
%
%   1D: [fi, fi_dx, fi_d2x] = binterp(f, xi, {borderCondition})
%   2D: [fi, fi_dx, fi_dy, fi_d2x, fi_d2y] = binterp(f, xi, yi, {borderCondition})
%
% Inputs: 
%                  f : input signal or image
%             xi, yi : interpolation coordinates (must be of same size)
%  {borderCondition} : 'mirror' (default) or 'periodic'
%
% Outputs: 
%                 fi : interpolated signal/image
%       fi_dx, fi_dy : 1st partial derivatives of input signal/image at interpolation coordinates
%     fi_d2x, fi_d2y : 2nd partial derivatives of input signal/image at interpolation coordinates
%
% For more information, see:
% [1] Unser, IEEE Signal Proc. Mag. 16(6), pp. 22-38, 1999
% [2] Unser et al., IEEE Trans. Signal Proc. 41(2), pp. 834-848, 1993

% Francois Aguet ? May 2012 (last modified May 11, 2012)
function [fi, fi_dx, fi_d2x] = binterp(f, xi, borderCondition) %#ok<STOUT,INUSD>