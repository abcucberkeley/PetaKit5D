%
% Examples for the binterp() function.
%
% Francois Aguet

%===============================================================================
% Perfect reconstruction for different border conditions
%===============================================================================
f = [4 6 8 10 9 1 2 8];
x = 1:numel(f);

figure;

f_rec = binterp(f, x, 'mirror');
subplot(2,1,1);
plot(x, f-f_rec);

f_rec = binterp(f, x, 'periodic');
subplot(2,1,2);
plot(x, f-f_rec);

%%
%===============================================================================
% 1D example with periodic boundary conditions
%===============================================================================

f = [0 3 2 1 0 3 2 1 0 3 2 1];
x = 1:numel(f);
xi = 1:0.01:numel(f);

[fi, fi_dx, fi_d2x] = binterp(f, xi, 'mirror');

figure;
hold on;
plot(x, f, 'ro');
plot(xi, fi, 'k-');
plot(xi, fi_dx, 'b-');
plot(xi, fi_d2x, 'c-');
legend('f[k]', 'f(x)', 'f''(x)', 'f''''(x)');
set(gca, 'XLim', [0.5 numel(f)+0.5]);

%% 
%===============================================================================
% Parametric curves
%===============================================================================
x_t = [1 0 -1 0];
y_t = [0 1 0 -1];

% closed curve
dt = 0.1;
t = 1:dt:length(x_t)+1-dt; % last node == first node

[fx, d_fx, d2_fx] = binterp(x_t, t, 'periodic');
[fy, d_fy, d2_fy] = binterp(y_t, t, 'periodic');

figure;
plot([x_t x_t(1)], [y_t y_t(1)], 'r.');
axis equal; hold on;
plot(fx, fy, 'ko');
plotcircle([0 0], 1, 'EdgeColor', 0.6*[1 1 1], 'LineStyle', '--');
axis(1.1*[-1 1 -1 1]);
legend('Samples/node points', 'Cubic B-spline', 'Circle', 'Location', 'NorthEastOutside');

%%
%===============================================================================
% 2D example with periodic boundary conditions
%===============================================================================
f = [0 3 2 1 0 3 2 1 0 3 2 1];
f = f'*f;
% f = f'*[0 3 2 1 0 3 2 1];
% f = f';
[ny,nx] = size(f);
[X,Y] = meshgrid(1:0.01:nx,1:0.01:ny);

[F, F_dx, F_dy, ~, ~] = binterp(f, X, Y, 'periodic');

figure;
subplot(2,2,1)
imagesc(f); colormap(gray(256)); axis image;
title('f');

subplot(2,2,2)
imagesc(F); colormap(gray(256)); axis image;
title('F');

subplot(2,2,3)
imagesc(F_dx); colormap(gray(256)); axis image;
title('d/dx F');

subplot(2,2,4)
imagesc(F_dy); colormap(gray(256)); axis image;
title('d/dy F');


xa = 1:0.01:nx;
ya = 1:0.01:ny;
[X,Y] = meshgrid(xa,ya);

tic;
F = binterp(f, X, Y, 'periodic');
t1 = toc;

[xm,ym] = meshgrid(1:nx,1:ny);
tic;
Fm = interp2(xm, ym, f, X, Y, 'spline');
t2 = toc;

figure;
subplot(1,2,1)
imagesc(F); colormap(gray(256)); axis image;
title(['F (' num2str(t1, '%.3f') ' s)']);

subplot(1,2,2)
imagesc(Fm); axis image; colormap(gray(256));
title(['Matlab spline (' num2str(t2, '%.3f') ' s)']);

%%

% test perfect reconstruction
nx = 20;
ny = 10;
S = rand(ny,nx);

xa = 1:nx;
ya = 1:ny;
[xi,yi] = meshgrid(xa,ya);
Srec = binterp(S, xi, yi, 'symmetric');
figure; imagesc(S-Srec); colormap(gray(256)); axis image; colorbar;

C = computeBSplineCoefficients(S);
k = [1 4 1]/6;
rec = conv2(k, k, padarrayXT(C, [1 1], 'symmetric'), 'valid');
figure; imagesc(rec-S); colormap(gray(256)); axis image; colorbar;




