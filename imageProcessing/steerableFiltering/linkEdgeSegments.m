% Francois Aguet, 10/2012

function [edgeMask] = linkEdgeSegments(edgeMask, linkIndex)

edgeMask = edgeMask~=0;
dims = size(edgeMask);
linkMask = zeros(dims);
for k = 1:size(linkIndex,1)
    [y0, x0] = ind2sub(dims, linkIndex(k,1));
    [y1, x1] = ind2sub(dims, linkIndex(k,2));
    iseg = bresenham([x0 y0], [x1 y1]);
    iseg = sub2ind(dims, iseg(:,2), iseg(:,1));
    linkMask(iseg) = 1;
end
linkMask = bwmorph(linkMask, 'thin');
edgeMask = edgeMask | linkMask;
