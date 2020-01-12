% function [inAngIdx] = angleFilter(vecx,vecy,ang_mid)
% filter in vectors between -pi/4 and pi/4 from vec_mid
function [inAngIdx] = angleFilter(vecx,vecy,vec_mid)
Npoints = length(vecx);
inAngIdx = false(Npoints,1);

for ii=1:Npoints
    vec = [vecx;vecy];
    ang = acos(vec'*vec_mid/(norm(vec)*norm(vec_mid)));
    if ang>-pi/3 && ang<pi/3
        inAngIdx(ii) = true;
    end
end
