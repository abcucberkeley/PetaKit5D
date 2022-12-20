function [imSize] = getImageSize(filePath)
% Get image size without loading the image
% 
% use Tiff for faster performance
% 
% Author: Xiongtao Ruan (01/2020)
% xruan (10/15/2020): add support for zarr format
% xruan (01/20/2021): use binary search to get the number of pages for tiff
% which is 10~20 times for large images faster than read each page linearly!
% xruan (01/20/2021): add option to use bioformat to get image size
% xruan (01/29/2021): check if jvm available, if not use binary search
% xruan (01/14/2022): use mex version by default


warning ('off','all');

if strcmp(filePath(end - 2 : end), 'tif') || strcmp(filePath(end - 3 : end), 'tiff')
    % tobj = Tiff(filePath, 'r');
    % Height = getTag(tobj, 'ImageLength');
    % Width = getTag(tobj, 'ImageWidth');
    
    % Zstack = 1;
%     while ~tobj.lastDirectory()
%         tobj.nextDirectory() ;
%         Zstack = Zstack + 1;
%     end
%    
    try
        imSize = getImageSize_mex(filePath);
    catch ME
        if ~usejava('jvm') && ~useParpool
            % use binary search 
            tobj = Tiff(filePath, 'r');     
            Height = getTag(tobj, 'ImageLength');
            Width = getTag(tobj, 'ImageWidth');

            s = 1;
            t = 1;
            is_valid = true;
            % first define (tight) upper bound of the stack number    
            while is_valid
                s = t;
                t = 16 * t;
                is_valid = set_tiff_directory_number(tobj, t); 
            end

            % narrow the range to find the stack number
            while s ~= t
                mid = round((s + t) / 2);
                is_valid = set_tiff_directory_number(tobj, mid);
                if is_valid
                    s = mid;
                else
                    t = mid - 1;
                end
            end

            Zstack = s;
        else
            % use bioFormat reader
            reader = bfGetReader(filePath); 
            Height = reader.getSizeY;
            Width = reader.getSizeX;
            Zstack = reader.getSizeT;
            if Zstack == 1
                obj = matlab.io.internal.BigImageTiffReader(filePath);
                Zstack = obj.NumImages;
            end

            reader.close()
        end
    %     close(tobj);
        imSize = [Height, Width, Zstack];
    end
elseif strcmp(filePath(end - 3 : end), 'zarr')
    try 
        bim = blockedImage(filePath, 'Adapter', CZarrAdapter);
    catch ME
        disp(ME);
        bim = blockedImage(filePath, 'Adapter', ZarrAdapter);
    end
    imSize = bim.Size;
end


end


function [is_valid] = set_tiff_directory_number(t, dirNum)

try 
    t.setDirectory(dirNum)
    is_valid = true;    
catch
    is_valid = false;
end

end


