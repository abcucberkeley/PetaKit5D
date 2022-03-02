function X=sepblockfun(X,blockdims,fun)
%Perform a separable operation on sub-blocks of an input array. Here, a function
%op() is said to be separable if for any array B with elements B(i,j,k,...), 
%the operation op(B(:)), resulting in a scalar, can be equivalently done by 
%applying op() first along i, then along j, then along k, etc...
%
%USAGE:
%
%      Y=sepblockfun(X,blockdims,fun)
%
%in:
%
%
%   X: A full array. If the ndSparse class defintion is on the path, then X 
%      can also be a regular sparse matrix  or ndSparse array. Performance might 
%      not be as strong as for full arrays, however.
%
%   blockdims: a vector of integers  specifying the dimensions of the
%              sub-blocks. The array X must partition evenly into blocks 
%              of this size. If blockdims(i) is set to Inf then it will be 
%              replaced with blockdims(i)=size(X,i).
%
%   fun:  function handle to an operation assumed to be separable
%         (Examples: max,min,sum,prod,mean, etc...). The function must
%         accept the input syntax fun(B,DIM) where B is an input array
%         and DIM is a dimension along which to operate.  Alternatively,
%         fun can be one of the following strings 'max','min','sum','mean',
%         'prod'.
%
%
%out:
%
%   Y: the output array. Y(i)=fun(Xi(:),1) where Xi is the i-th sub-block of
%      the input array X.
%
%
%EXAMPLE 1: Divide a 400x400x400 array into 10x10x10 blocks. Return the blockwise, 
%mean, max, and min. of each block, each organized as a 40x40x40 array.
%
%   A=rand(400,400,400);
%   Ameans=sepblockfun(A,[10,10,10],@mean);
%   Amins=sepblockfun(A,[10,10,10],'min' ); 
%   Amaxs=sepblockfun(A,[10,10,10], @(B,d) max(B,[],d)  );
%
%EXAMPLE 2: Not all operations satisfy the separability property, but
%sometimes inseparable operations can be decomposed into separable ones. As
%an example, we take the blockwise standard deviations of the  same array
%from Example 1.
%
%   Astds=sqrt( sepblockfun(A.^2,[10,10,10],'mean') - Ameans.^2 );  
%
%
% Written by Matt Jacobson, 2014

if issparse(X)&& exist('ndSparse','class') && ~isa(X,'ndSparse')

     X=ndSparse(X);
    
end


if ischar(fun)

  switch fun

    case 'max'
     
         fun=@(b,d) max(b,[],d);

    case 'min'

         fun=@(b,d) min(b,[],d);

    case 'sum'

         fun=@sum;

    case 'mean'

         fun=@mean;

    case 'prod'

        fun=@prod;


    otherwise

     error 'Unrecognized fun() selection'

  end


end


nn=max(length(blockdims),ndims(X));
blockdims(end+1:nn)=1;

[sz{1:nn}]=size(X); %M is the original array
sz=[sz{:}];

idx=~isfinite(blockdims);
blockdims(idx)=sz(idx);

newdims=sz./blockdims;

args=num2cell([blockdims;newdims]);

X=reshape(X,args{:});

for ii=1:nn


 X=fun(X,2*ii-1);
 
end

X=reshape(X,newdims);
