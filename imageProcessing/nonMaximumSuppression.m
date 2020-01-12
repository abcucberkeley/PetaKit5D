%nms = nonMaximumSuppression(res, th)
%
% Inputs:   res : response
%            th : orientation
%
% Uses the grid conventions of steerableDetector()

% Francois Aguet

function res = nonMaximumSuppression(res, th)

[ny,nx] = size(res);

res = padarrayXT(res, [1 1], 'symmetric');

[x,y] = meshgrid(1:nx,1:ny);

% +1 interp
A1 = interp2(res, x+1+cos(th), y+1+sin(th),'linear',0);

% -1 interp
A2 = interp2(res, x+1-cos(th), y+1-sin(th),'linear',0);

res = res(2:end-1,2:end-1);

res(res<A1 | res<A2) = 0;
