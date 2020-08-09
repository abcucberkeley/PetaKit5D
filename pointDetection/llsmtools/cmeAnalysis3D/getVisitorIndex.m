%[vidx] = getVisitorIndex(lftData, mCh) identifies trajectories corresponding to objects 'visiting' the TIRF field

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
