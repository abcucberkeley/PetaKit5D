classdef CZarrAdapter < images.blocked.Adapter
    % C-zarr based zarr adapter

    % Author: Xiongtao Ruan (07/25/2022)


    properties (SetAccess = protected, GetAccess = public)
        ZarrObj
        ZarrInfo
    end
     
    % Read methods
    methods
        function openToRead(obj, loc)
            if ~isfile([loc + "/.zgroup"])
                obj.ZarrObj.loc = loc;
            else
                if exist(loc + "/L_1_1_1", 'dir')
                    obj.ZarrObj.loc = loc + "/L_1_1_1";
                else
                    error('Unable to open %s', loc);
                end
            end
            
            % parse the .zarray file
            zarrayFn = [obj.ZarrObj.loc + '/.zarray'];
            fid = fopen(zarrayFn);
            raw = fread(fid);
            fclose(fid);
            str = char(raw');
            jdata = jsondecode(str);
            
            % assign val to obj
            obj.ZarrObj.chunks = jdata.chunks(:)';
            obj.ZarrObj.shape = jdata.shape(:)';
            if contains(jdata.dtype, 'u2')
                obj.ZarrObj.dtype = 'uint16';
            elseif contains(jdata.dtype, 'f4')
                obj.ZarrObj.dtype = 'single';
            end        
        end
         
        function info = getInfo(obj)
            info.Size = obj.ZarrObj.shape;
            info.IOBlockSize = obj.ZarrObj.chunks;
            % xruan (11/18/2020): add support for multiple data types
            switch string(obj.ZarrObj.dtype)
                case {"float32", "single"}
                    info.Datatype = "single";
                case "uint16"
                    info.Datatype = "uint16";
            end
            % info.Datatype = string(obj.ZarrObj.dtype.name);
            info.InitialValue = cast(0, info.Datatype);
            info.UserData = [];
            % For use later in getIOBlock
            obj.ZarrInfo =  info;
        end
         
        function data = getIOBlock(obj, ioBlockSub, level)
            assert(level==1)
            regionStart = (ioBlockSub-1).*obj.ZarrInfo.IOBlockSize + 1;
            regionEnd = min(obj.ZarrObj.shape, regionStart + obj.ZarrInfo.IOBlockSize - 1);
             
            data = readzarr(obj.ZarrObj.loc, 'inputBbox', [regionStart, regionEnd]);
        end
         
        function data = getIORegion(obj, regionStart, regionEnd) 
            bbox = [regionStart(:)', regionEnd(:)'];
            data = readzarr(obj.ZarrObj.loc, 'inputBbox', bbox);
        end
    end
     
    % Write methods
    methods
        function openToWrite(obj, loc, info, level)
            assert(level==1)
            if numel(info.Size) == 2
                % info.Size = [info.Size, 1];
                % info.IOBlockSize = [info.IOBlockSize, 1];
            end
            
            if ispc && numel(char(loc)) > 200
                [pth, loc] = fileparts(loc);
                cd(pth);                    
            end

            dimSeparator = '.';
            if prod(ceil(info.Size / info.IOBlockSize)) > 10000
                dimSeparator = '/';
            end

            createzarr(char(loc), dataSize=info.Size, blockSize=info.IOBlockSize, ...
                dtype=info.Datatype, order='F', compressor='zstd', zarrSubSize=[], ...
                dimSeparator=dimSeparator);

            openToRead(obj, loc);
            obj.ZarrInfo =  info;
        end
        

        function openInParallelToAppend(obj, loc)
            obj.ZarrObj = openToRead(obj, loc);
            obj.ZarrInfo = obj.getInfo();
        end

        function setIOBlock(obj, ioBlockSub, level, data)
            assert(level==1)
            if numel(ioBlockSub) == 2
                ioBlockSub = [ioBlockSub, 1];
            end
            regionStart = (ioBlockSub-1).*obj.ZarrInfo.IOBlockSize + 1;
            regionEnd = min(obj.ZarrObj.shape, regionStart + obj.ZarrInfo.IOBlockSize - 1);
            % say a 4x1 data as 2D in python. Error looks like:
            writezarr(data, obj.ZarrObj.loc, 'blockSize', obj.ZarrObj.chunks, ...
                'bbox', [regionStart, regionEnd], create=false);
        end
        
        function setRegion(obj, regionStart, regionEnd, data)
            % TODO - handle singleton dimensions, not sure how to convert
            writezarr(data, obj.ZarrObj.loc, 'blockSize', obj.ZarrObj.chunks, ...
                'bbox', [regionStart, regionEnd], create=false);
        end
        
        function setData(obj, data)
            % write the whole volume with multi-processing
            % TODO - handle singleton dimensions, not sure how to convert
            writezarr(data, obj.ZarrObj.loc, 'blockSize', obj.ZarrObj.chunks, ...
                'bbox', [1, 1, 1, size(data)]);
        end        
    end
end
