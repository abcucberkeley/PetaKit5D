function nHood = bwnHood3D(conn)
%NHOOD3D creates basic 6,18 and 26 connectivity 3x3x3 3D neighborhoods for morphological operations
%
% nHood = nHood3D;
%
% nHood = nHood3D(conn);
% 
% Input: 
% 
%   conn - Either 6, 18 or 26. This number specifies the type of
%   neighborhood to return. Optional. Default is 26.
% 
% 
% Output:
%
%   nHood - A 3x3x3 binary matrix with 6,18 or 16 1s as specified by 'conn'
%
%
%
% Hunter Elliott
% 4/2010
%

if nargin < 1 || isempty(conn)    
    conn = 26;    
end

nHood = false([3 3 3]);

switch conn
    
    case 6        
        
        nHood([5 11 13 15 17 23]) = true;
        
    case 18
        
        nHood([5 10:18 23]) = true;
        
    case 26
        
        nHood([1:13 15:end]) = true;
        
    otherwise
        nHood = [];
        error('Invalid input connectivity number: must be 6, 18 or 26')
        
end


        