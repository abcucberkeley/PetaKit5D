function thinner = bwthin(bw)
%BWTHIN thins a 2D or 3D binary image
% 
% thinner = bwthin(bw)
% 
% This thins the input binary 2D or 3D image with the use of a hit-miss
% operation and the standard structuring elements for binary thinning.
% 
% NOTE: This function is slow! Should really only be used for 3D if needed.
% For 2D thinning, use bwmorph instead!
%
% Hunter Elliott
% 3/2010


%Create the structuring elements for the hit miss operation in 2D and 3D.
if ndims(bw) == 2
    %Neighborhoods
    n1   = [-1 -1 -1 ;...
             0  1  0 ; ...
             1  1  1];

    n2    = [ 0 -1 -1 ;...
              1  1 -1;...
              0  1  0];
          
    nOuterLoop = 1;
elseif ndims(bw) == 3
    
   n1(:,:,1) = [ -1 -1 -1;...
                  0  0  0;...
                  1  1  1];
            
              
   n1(:,:,2) = [-1 -1 -1;...
                 0  1  0;...
                 1  1  1];
            
   n1(:,:,3) = [-1 -1 -1;...
                 0  0  0;...
                 1  1  1];
        
        
        
   n2(:,:,1) = [  0 -1 -1;...
                  1  0 -1;...
                  0  1  0];
            
   n2(:,:,2) = [  0 -1 -1;...
                  1  1 -1;...
                  0  1  0];
            
   n2(:,:,3) = [  0 -1 -1;...
                  1  0 -1;...
                  0  1  0];  
              
  nOuterLoop = 4;
        
end       
      
      
thinner = bw;      

%This is certainly not the fastest way to do this... but it works...

%Apply the hit-miss operation using the neighborhoods and their 90degree
%rotations.
for l = 1:nOuterLoop

    for j = 1:4

        thinner = ~bwhitmiss(thinner,n1) & thinner;

        thinner = ~bwhitmiss(thinner,n2) & thinner;

        for k = 1:size(n1,3)
            n1(:,:,k) = rot90(n1(:,:,k));
            n2(:,:,k) = rot90(n2(:,:,k));
        end

    end

    for k = 1:3
        n1(k,:,:) = rot90(squeeze(n1(k,:,:)));
        n2(k,:,:) = rot90(squeeze(n2(k,:,:)));
    end
end
