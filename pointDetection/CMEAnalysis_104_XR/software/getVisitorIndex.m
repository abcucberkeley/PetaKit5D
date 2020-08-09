%[vidx] = getVisitorIndex(lftData, mCh) identifies trajectories corresponding to objects 'visiting' the TIRF field
%
% Copyright (C) 2017, Danuser Lab - UTSouthwestern 
%
% This file is part of CMEAnalysis_Package.
% 
% CMEAnalysis_Package is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% CMEAnalysis_Package is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with CMEAnalysis_Package.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

function vidx = getVisitorIndex(lftData, mCh)

if nargin<2
    mCh = 1;
end

minLength = min(arrayfun(@(i) size(i.gapMat_Ia,2), lftData));

A = arrayfun(@(i) i.A(:,1:minLength,mCh), lftData, 'UniformOutput', false);
A = vertcat(A{:});
lftV = arrayfun(@(i) i.lifetime_s(i.catIdx==1), lftData, 'unif', 0);
lifetime_s_all = vertcat(lftV{:});

% 95th percentile of the reference (above threshold) distribution
tx = 30; % not sensitive, tested on 1 and 0.5 frame/sec data
pRef = prctile(A(lifetime_s_all>=tx,1:3,mCh), 95, 1);

nd = numel(lftData);
vidx = cell(1,nd);
for i = 1:nd
    vidx{i} = sum(lftData(i).A(:,1:3,mCh)>repmat(pRef, [size(lftData(i).A,1) 1]),2)>0 & lftV{i}<tx;
end
