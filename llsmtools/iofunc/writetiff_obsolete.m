%writetiff(img, filepath, varargin) writes a TIFF stack using libtiff
% Stores TIFFs as 64-bit

% Francois Aguet, 05/21/2013

function writetiff(img, filepath, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('Compression', 'lzw', @(x) any(strcmpi(x, {'none', 'lzw'})));
ip.addParamValue('Mode', 'libtiff', @(x) any(strcmpi(x, {'libtiff', 'imwrite'})));
ip.parse(varargin{:});

t = Tiff(filepath,'w');

[ny,nx,N] = size(img);
w = whos('img');

switch lower(ip.Results.Mode)
    case 'libtiff'
        tagstruct.ImageLength = ny;
        tagstruct.ImageWidth = nx;
        tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
        switch lower(ip.Results.Compression)
            case 'none'
                tagstruct.Compression = Tiff.Compression.None;
            case 'lzw'
                tagstruct.Compression = Tiff.Compression.LZW;
        end
        tagstruct.BitsPerSample = 8*w.bytes/nx/ny/N;
        tagstruct.SamplesPerPixel = 1;
        tagstruct.RowsPerStrip = 16;
        if any(strcmpi(w.class, {'uint8', 'uint16', 'uint32', 'uint64'}))
            tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
        elseif any(strcmpi(w.class, {'int8', 'int16', 'int32', 'int64'}))
            tagstruct.SampleFormat = Tiff.SampleFormat.Int;
        else
            tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
        end
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tagstruct.SubFileType = Tiff.SubFileType.Page;
        tagstruct.PageNumber = [1 N];
        tagstruct.Software = 'MATLAB';
        
        for k = 1:N;
            t.setTag(tagstruct);
            t.write(img(:,:,k));
            t.writeDirectory();
        end
        t.close();
    case 'imwrite'
        imwrite(img(:,:,1), filepath, 'tif', 'compression' , ip.Results.Compression);
        for k = 2:N
            imwrite(img(:,:,k), filepath, 'tif', 'compression' , ip.Results.Compression, 'writemode', 'append');
        end
end
