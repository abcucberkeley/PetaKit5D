function B = tls(X,Y,useGPU)
%Naive total least squares taken from wikipedia
%[m, n]   = size(X);            % n is the width of X (X is m by n)
[~, n]   = size(X);            % n is the width of X (X is m by n)
Z       = [X Y];              % Z is X augmented with Y.
if(useGPU)
    Z = gather(Z);          
    [~, ~, V] = svd(Z,0);           % find the SVD of Z.
    V = gpuArray(V);
else 
    [~, ~, V] = svd(Z,0); 
end
VXY     = V(1:n,1+n:end);     % Take the block of V consisting of the first n rows and the n+1 to last column
VYY     = V(1+n:end,1+n:end); % Take the bottom-right block of V.
B       = -VXY/VYY;
end