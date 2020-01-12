%[stack] = insertMissingFrame(stack, frame) uses normalized cross-correlation to identify the position of the dropped frame in the stack

% Francois Aguet, 10/30/2013

function [stack] = insertMissingFrame(stack, frame)

nz = size(stack,3);
cv = zeros(1,nz);
for i = 1:nz
    tmp = stack(:,:,i);
    % norm. corr:
    cv(i) = sum(frame(:).*tmp(:))/sqrt(sum(tmp(:).^2)*sum(frame(:).^2));
end
[~,maxIdx] = max(cv);

% figure; plot(cv);

if maxIdx==1
    stack = cat(3, frame, stack);
elseif maxIdx==nz
    stack = cat(3, stack, frame);
elseif cv(maxIdx-1)>cv(maxIdx+1)
    stack = cat(3, stack(:,:,1:maxIdx-1), frame, stack(:,:,maxIdx:end));
else
    stack = cat(3, stack(:,:,1:maxIdx), frame, stack(:,:,maxIdx+1:end));
end
