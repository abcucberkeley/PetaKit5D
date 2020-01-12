function [dY,dX]=createDiffMatrix(Pi,Pg)
% createDiffMatrix is an accessory C-MEX function for vectorFieldDiv
%
% SYNOPSIS   [dX,dY]=createDiffMatrix(Pi,Pg)
%
% INPUT      Pi and Pg are the matrices containing the set of 2D point coordinates.
%
%                   M=[ y1 x1     and   N=[ y1 x1
%                       y2 x2              y2 x2
%                        ...                ...
%                       ym xm ]            yn xn ]
%
% OUTPUT   dY :
%          dX :
% 
% REMARK   
%
% C-MEX file - Aaron Ponti 11/26/02
