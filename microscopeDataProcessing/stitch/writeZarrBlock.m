function [] = writeZarrBlock(bim, blockSub, data, level, Mode) 
% write block in write or read mode

if strcmp(Mode, 'w')
    bim.setBlock(blockSub, data);
else
    % adapted from blockedImage.setBlock
    if any(blockSub == bim.SizeInBlocks)
        % Edge block, check if data needs to be trimmed to fit 
        totalData = (blockSub-1).*bim.BlockSize + size(data, 1:bim.NumDimensions);
        trimAmount = max(0, totalData - bim.Size);
        if any(trimAmount)
            indStruct.type = '()';
            indStruct.subs = cell(1,bim.NumDimensions);
            for dInd = 1:bim.NumDimensions
                indStruct.subs{dInd} = 1:(size(data,dInd)-trimAmount(dInd));
            end
            data = subsref(data, indStruct);
        end
    end

    bim.Adapter.setIOBlock(blockSub, level, data);
end


end
