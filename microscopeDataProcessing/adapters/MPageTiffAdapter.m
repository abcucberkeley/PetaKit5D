classdef MPageTiffAdapter < images.blocked.Adapter
% adapter for 3D tiff reading (multipage tiff)
%   
% adapted from Tiff adapter and AOLLSMAdater
% 
% Author: Xiongtao Ruan (10/02/2020)
% xruan (10/10/2020): add support for writing
% xruan (12/23/2021): add support for 2d image
% xruan (07/25/2022): change to use parallelWriteTiff as default write method

    properties
        Extension (1,1) string = "tiff"
        Compression = Tiff.Compression.LZW
    end
    
    properties (Access = private)
        FileName (1,1) string
        ITRObject = [];   % For read
        TIFFObject = [];  % For write
        ActiveLevel
        SizeInBlocks
        Info
        TIFFInfoCache = containers.Map()
        LastWrittenLevel = 1
    end    

    % Read methods
    methods
        function openToRead(obj, tiffFileName)
            obj.FileName = tiffFileName;
            
            try
                obj.ITRObject = matlab.io.internal.BigImageTiffReader(obj.FileName);
            catch ME
                newException = MException('images:bigimages:couldNotReadUserData', ME.message);
                throw(newException)
            end
            
            imageSize = [obj.ITRObject.ImageHeight, obj.ITRObject.ImageWidth, obj.ITRObject.NumImages];
            ioBlockSize = [obj.ITRObject.ImageHeight, obj.ITRObject.ImageWidth, 1];
            datatype = obj.ITRObject.MLType;
            
            obj.Info.Size = imageSize;
            obj.Info.IOBlockSize = ioBlockSize;
            obj.Info.Datatype = datatype;            
            obj.Info.InitialValue = cast(0, obj.Info.Datatype(1));                       
            % obj.Info.UserData = images.blocked.internal.loadUserData(obj.FileName);
            
            obj.SizeInBlocks = ceil(imageSize./ioBlockSize);
            
            % Find the 'finest' level (most pixels)
            numPixels = prod(imageSize,3);
            if imageSize(3) == 1
                finestLevel = 1;
            else
                [~, finestLevel] = max(numPixels);
            end
            
            % Re-open the finest level by default
            obj.ITRObject = matlab.io.internal.BigImageTiffReader(...
                obj.FileName,"ImageIndex", finestLevel);
            obj.ActiveLevel = finestLevel;
        end
        
        function info = getInfo(obj)
            info = obj.Info;
        end
        
        function data = getIOBlock(obj, ioBlockSub, level)
            assert(level==1, "Only single level data is supported");
            
            % zInd = (ioBlockSub(3)-1)*obj.Info.IOBlockSize(3) + 1 : ioBlockSub(3)*obj.Info.IOBlockSize(3);
            zInd = ioBlockSub(3);
            tiffFileName = obj.FileName;
            
            %--------------------------------------------------------------
            if ~obj.TIFFInfoCache.isKey(tiffFileName)                
                % We have not yet read the info struct for this file, read
                % and cache. 
                %----------------------------------------------------------
                %This cache can get large, its worth exploring
                % if we can skip imread/imfinfo etc and just plain fread
                % the data (this will only work for uncompressed TIFF files
                % and might be significantly faster).
                tinfo = imfinfo(tiffFileName);
                obj.TIFFInfoCache(tiffFileName) = tinfo;
                %----------------------------------------------------------
            end

            try
                % data = readtiff(tiffFileName, "range", zInd, "info", obj.TIFFInfoCache(tiffFileName));
                data = imread(tiffFileName, "Index", zInd, ...
                    "Info", obj.TIFFInfoCache(tiffFileName));
                % data = data(yInd, xInd, :);
            catch
                % In case the data is corrupt
                warning("Could not read slice: " + num2str(zInd) + " from file " + tiffFileName);
                data = zeros(obj.Info.IOBlockSize,'uint16');
            end
        end

        function data = getIORegion(obj, regionStart, regionEnd)
            tiffFileName = obj.FileName;            
            data = readtiff(char(tiffFileName), [regionStart(3), regionEnd(3)]);
            try
                data = crop3d_mex(data, [regionStart(1 : 2), 1, regionEnd(1 : 2), size(data, 3)]);
            catch ME
                disp(ME);
                data = data(regionStart(1) : regionEnd(1), regionStart(2) : regionEnd(2), :);
            end
        end        
    end        
        
    methods
        function  openToWrite(obj, destination, info, level)
            % Write mode, always write bigTiff, since we are usually
            % dealing with large data. And always interleaved (chunky)
            % since thats the Tiff baseline spec.
            % assert(size(info.Size, 2)==2 || size(info.Size, 2)==3 ,...
            %     "Only Grayscale and RGB TIFF files are supported");
            
            if level<obj.LastWrittenLevel
                error(message('images:blockedImage:TIFFLevelIsClosed'))                
            end
            obj.LastWrittenLevel = level;
            
            obj.Info = info;            
            obj.FileName = destination;
            try
                if isfile(destination)
                    assert(~isempty(obj.TIFFObject), "File exists, but not writing multiple levels");
                    obj.TIFFObject.writeDirectory();
                else
                    path = fileparts(destination);
                    images.blocked.internal.createFolder(path);
                    obj.TIFFObject = Tiff(destination, 'w');
                end
            catch WRITEEXP
                EXP = MException('images:bigimage:couldNotCreate',...
                    message('images:bigimage:couldNotCreate', destination));
                EXP = EXP.addCause(WRITEEXP);
                throw(EXP)
            end
            
            % images.blocked.internal.saveUserData(obj.FileName, info.UserData);
            
            % Setup the tags
            tagStruct.ImageLength = info.Size(level,1);
            tagStruct.ImageWidth = info.Size(level,2);
            tagStruct.RowsPerStrip = 512; % http://www.awaresystems.be/imaging/tiff/tifftags/rowsperstrip.html

            if size(info.Size,2)==2 || info.Size(level,3) == 1
                tagStruct.Photometric = Tiff.Photometric.MinIsBlack;
                tagStruct.SamplesPerPixel = 1;
                tagStruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            elseif info.Size(level,3) == 3
                tagStruct.Photometric = Tiff.Photometric.RGB;
                tagStruct.SamplesPerPixel = 3;
                tagStruct.PlanarConfiguration = Tiff.PlanarConfiguration.Separate;
            elseif info.Size(level,3) > 3
                tagStruct.Photometric = Tiff.Photometric.MinIsBlack;                
                tagStruct.SamplesPerPixel = 1;
                tagStruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            else
                error(message('images:bigimage:tiffChannelCount'))
            end
            
            switch info.Datatype(level)
                case 'logical'
                    tagStruct.BitsPerSample = 1;
                case {'uint8', 'int8'}
                    tagStruct.BitsPerSample = 8;
                case {'uint16', 'int16'}
                    tagStruct.BitsPerSample = 16;
                case {'uint32', 'int32', 'single'}
                    tagStruct.BitsPerSample = 32;
                case {'double'}
                    tagStruct.BitsPerSample = 64;
                otherwise
                    assert(false)
            end
            
            switch info.Datatype(level)
                case {'logical', 'uint8', 'uint16', 'uint32'}
                    tagStruct.SampleFormat = Tiff.SampleFormat.UInt;
                case {'int8', 'int16', 'int32'}
                    tagStruct.SampleFormat = Tiff.SampleFormat.Int;
                case {'single','double'}
                    tagStruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
                otherwise
                    assert(false)
            end
            
            % for big data, save as bigtiff
            if prod(info.Size) * tagStruct.BitsPerSample / 8 > 2^32-1
                obj.TIFFObject = Tiff(destination, 'w8');
            end
            
            % tagStruct.TileWidth = info.IOBlockSize(level,2);
            % tagStruct.TileLength = info.IOBlockSize(level,1);
            tagStruct.Photometric = obj.TIFFObject.Photometric.MinIsBlack;
            tagStruct.Compression = obj.Compression;
            tagStruct.Software = 'MATLAB blockedImage';
            
            obj.Info.tagStruct = tagStruct;
            obj.SizeInBlocks = ceil(info.Size./info.IOBlockSize);
        end
        
        function setIOBlock(obj, ioBlockSub, level, data)
            assert(level==1)            
            zRange = (ioBlockSub(3)-1).*obj.Info.IOBlockSize(level,3) + (1 : size(data, 3)); 
            % data = permute(data, [1, 2, 4, 3]);
            for zInd = 1 : size(data, 3)
                obj.TIFFObject.setTag(obj.Info.tagStruct);                      
                obj.TIFFObject.write(data(:, :, zInd));
                if zRange(zInd) ~= obj.Info.Size(3)
                   obj.TIFFObject.writeDirectory();
                end
            end
        end

        function setRegion(obj, regionStart, regionEnd, data)
            disp('Only support the write of whole data for tiff')
            setData(obj, data);
        end
        
        function setData(obj, data)
            writetiff(data, obj.FileName);            
        end        
        
        function close(obj)
            if ~isempty(obj.TIFFObject)
                % Close the writer object
                obj.TIFFObject.close();
            end
        end
    end
end
