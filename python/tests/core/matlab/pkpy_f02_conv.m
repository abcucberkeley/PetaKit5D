function y = pkpy_f02_conv(x, cls)
% class functions called directly, as PetaKit5D's code calls them (not feval or cast)
switch cls
    case 'int8'
        y = int8(x);
    case 'uint8'
        y = uint8(x);
    case 'int16'
        y = int16(x);
    case 'uint16'
        y = uint16(x);
    case 'int32'
        y = int32(x);
    case 'uint32'
        y = uint32(x);
    case 'int64'
        y = int64(x);
    case 'uint64'
        y = uint64(x);
    case 'single'
        y = single(x);
    case 'double'
        y = double(x);
    case 'logical'
        y = logical(x);
end
end
