classdef ZarrAdapter < images.blocked.Adapter
    %

    %   Copyright 2020 The MathWorks, Inc.
    % 
    % xruan (11/12/2021): add support for hierachical group zarr file, for
    % now, we read the 'L_1_1_1' one if it is a group
    % xruan (11/14/2021): add methods to directly read/write regions for
    % given start and end coordinates, which is much faster than reading
    % block by block
    % xruan (03/11/2022): change to F order for zarr file
    
    properties (SetAccess = protected, GetAccess = public)
        ZarrObj
        ZarrInfo
    end
     
    % Read methods
    methods
        function openToRead(obj, loc)
            if ~isfile([loc + "/.zgroup"])
                obj.ZarrObj = py.zarr.open(loc);
            else
                if exist(loc + "/L_1_1_1", 'dir')
                    obj.ZarrObj = py.zarr.open(loc + "/L_1_1_1");
                else
                    error('Unable to open %s', loc);
                end
            end
        end
         
        function info = getInfo(obj)
            info.Size = cellfun(@double, cell(obj.ZarrObj.shape));
            info.IOBlockSize = cellfun(@double, cell(obj.ZarrObj.chunks));
            % xruan (11/18/2020): add support for multiple data types
            switch string(obj.ZarrObj.dtype.name)
                case "float32"
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
            % In python indices
            regionStart = uint32((ioBlockSub-1).*obj.ZarrInfo.IOBlockSize);
            regionEnd = uint32(double(regionStart) + obj.ZarrInfo.IOBlockSize);
            pydata = py.zarrAPI.getZarrRegion(obj.ZarrObj, regionStart, regionEnd);
             
            data = cast(pydata, obj.ZarrInfo.Datatype);
        end
         
        function data = getIORegion(obj, regionStart, regionEnd)
            % In python indices
            regionStart = uint32(regionStart - 1);
            regionEnd = uint32(regionEnd);
            pydata = py.zarrAPI.getZarrRegion(obj.ZarrObj, regionStart, regionEnd);
             
            data = cast(pydata, obj.ZarrInfo.Datatype);
        end
        
    end
     
    % Write methods
    methods
        function openToWrite(obj, loc, info, level)
            assert(level==1)
            % TODO - map info.Datatype to corresponding Zarr type.
            % xruan (11/18/2020): add support for multiple data types
            switch info.Datatype
                case 'single'
                    dtype = 'f4';
                case 'uint16'
                    dtype = 'u2';
            end
            if numel(info.Size) == 2
                % info.Size = [info.Size, 1];
                % info.IOBlockSize = [info.IOBlockSize, 1];
            end
            
            if ispc && numel(char(loc)) > 200
                [pth, loc] = fileparts(loc);
                cd(pth);                    
            end
            
            obj.ZarrObj = py.zarr.open(loc,...
                pyargs('mode', 'w',...
                'shape', uint32(info.Size),...
                'chunks', uint32(info.IOBlockSize),...
                'order', 'F',...
                'dtype', dtype));
            obj.ZarrInfo =  info;
        end
        
%         function openInParallelToAppend(obj, loc)
%             % assert(level==1)
%             % TODO - map info.Datatype to corresponding Zarr type.
% %             obj.ZarrObj = py.zarr.open(loc,...
% %                 pyargs('mode', 'w',...
% %                 'shape', uint32(info.Size),...
% %                 'chunks', uint32(info.IOBlockSize),...
% %                 'dtype', 'u2')); % uint16 hardcoded
%             obj.ZarrObj = py.zarr.open(loc);
%             % obj.ZarrInfo =  info;
%             obj.ZarrInfo = obj.getInfo();
% 
%         end
        function openInParallelToAppend(obj, loc)
            obj.ZarrObj = py.zarr.open(loc);
            obj.ZarrInfo = obj.getInfo();
        end

        function setIOBlock(obj, ioBlockSub, level, data)
            assert(level==1)
            % In python indices
            if numel(ioBlockSub) == 2
                ioBlockSub = [ioBlockSub, 1];
            end
            regionStart = uint32((ioBlockSub-1).*obj.ZarrInfo.IOBlockSize);
            regionEnd = uint32(double(regionStart) + obj.ZarrInfo.IOBlockSize);
            % TODO - handle singleton dimensions, not sure how to convert
            % say a 4x1 data as 2D in python. Error looks like:
            %    Python Error: ValueError: parameter 'value': expected array with shape (4,1), got (4,)
            pydata = py.numpy.array(data);
            py.zarrAPI.setZarrRegion(obj.ZarrObj, regionStart, regionEnd, pydata);
        end
        
        function setRegion(obj, regionStart, regionEnd, data)
            % In python indices
            regionStart = uint32(regionStart-1);
            regionEnd = uint32(regionEnd);
            % TODO - handle singleton dimensions, not sure how to convert
            % say a 4x1 data as 2D in python. Error looks like:
            %    Python Error: ValueError: parameter 'value': expected array with shape (4,1), got (4,)
            pydata = py.numpy.array(data);
            py.zarrAPI.setZarrRegion(obj.ZarrObj, regionStart, regionEnd, pydata);
        end
        
        function setData(obj, data)
            % write the whole volume with multi-processing
            % In python indices
            % TODO - handle singleton dimensions, not sure how to convert
            % say a 4x1 data as 2D in python. Error looks like:
            %    Python Error: ValueError: parameter 'value': expected array with shape (4,1), got (4,)
            pydata = py.numpy.array(data);
            py.zarrAPI.setZarrData(obj.ZarrObj, pydata);            
        end
    end
end
