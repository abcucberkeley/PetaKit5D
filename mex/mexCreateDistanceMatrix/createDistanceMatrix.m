function D=createDistanceMatrix(M,N)
% createDistanceMatrix calculates the distance matrix for two sets of points
%
% SYNOPSIS   D=createDistanceMatrix(M,N)
%
% INPUT      M and N are the matrices containing the set of point coordinates.
%            M and N can represent point positions in 1, 2 and 3D, as follows.
%            
%            In 1D: M=[ x1        and   N=[ x1
%                       x2                  x2
%                       ...                ... 
%                       xm ]                xn ]
%
%            In 2D:
%                   M=[ y1 x1     and   N=[ y1 x1
%                       y2 x2              y2 x2
%                        ...                ...
%                       ym xm ]            yn xn ]
%
%            In 3D:
%                   M=[ y1 x1 z1  and   N=[ y1 x1 z1
%                       y2 x2 z2            y2 x2 z2
%                         ...                ...
%                       ym xm zm ]          yn xn zn ]
%
%
% OUTPUT   D : distance matrix D=(dij), i=1..m, j=1..n
% 
% REMARK   For 1D, both positive and negative distances are returned.
%
% C-MEX file - Aaron Ponti 28/08/15
