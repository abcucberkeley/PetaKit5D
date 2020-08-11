function [goodNum] = findGoodFactorNumber(givenNum, directionFlag, allowOddSeven)
% find the nearest number that can only be factorized by 2, 3, 5 or 7. 
% 
% Xiongtao Ruan (03/10/2020)
% xruan (03/12/2020): in some cases, cufft fails even for good numbers. It
% seems an even number with factor 7 works.
% xruan (06/25/2020): add option for all odd number with factor 7.

if nargin < 3
    allowOddSeven = false;
end

if nargin < 2
    directionFlag = 1;
end

% givenFactors = [2, 3, 5, 7]';
goodNum = givenNum;

primeFactors = factor(goodNum);
while primeFactors(end) > 7 || (~allowOddSeven && (primeFactors(1) > 2 && primeFactors(end) == 7))
    goodNum = goodNum + directionFlag;
    primeFactors = factor(goodNum);
end

end
