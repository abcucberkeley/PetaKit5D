function [byteNum] = dataTypeToByteNumber(dtype)
% get byte number for the given data type in matlab


switch dtype
    case 'uint8'
        byteNum = 1;
    case 'uint16'
        byteNum = 2;
    case 'single'
        byteNum = 4;
    case 'double'
        byteNum = 8;
    otherwise
        error('Unsupported data type %s', dtype);
end

end