classdef N5Adapter < images.blocked.Adapter
    % N5 adapter
    % adapted from ZarrAdapter.m

    % Author: Xiongtao Ruan (07/28/2023)
    
    
    properties (SetAccess = protected, GetAccess = public)
        N5Obj
        N5Info
    end
     
    % Read methods
    methods
        function openToRead(obj, loc)
            if ~isfile([loc + "/.zgroup"])
                store = py.zarr.N5Store(loc);
                obj.N5Obj = py.zarr.open(store);
            else
                if exist(loc + "/s0", 'dir')
                    store = py.zarr.N5Store(loc + "/s0");
                    obj.N5Obj = py.zarr.open(store);
                else
                    error('Unable to open %s', loc);
                end
            end
        end
         
        function info = getInfo(obj)
            info.Size = cellfun(@double, cell(obj.N5Obj.shape));
            info.IOBlockSize = cellfun(@double, cell(obj.N5Obj.chunks));
            % xruan (11/18/2020): add support for multiple data types
            switch string(obj.N5Obj.dtype.name)
                case "float32"
                    info.Datatype = "single";
                case "uint16"
                    info.Datatype = "uint16";
            end
            % info.Datatype = string(obj.N5Obj.dtype.name);
            info.InitialValue = cast(0, info.Datatype);
            info.UserData = [];
            % For use later in getIOBlock
            obj.N5Info =  info;
        end
         
        function data = getIOBlock(obj, ioBlockSub, level)
            assert(level==1)
            % In python indices
            regionStart = uint32((ioBlockSub-1).*obj.N5Info.IOBlockSize);
            regionEnd = uint32(double(regionStart) + obj.N5Info.IOBlockSize);
            pydata = py.zarrAPI.getZarrRegion(obj.N5Obj, regionStart, regionEnd);
             
            data = cast(pydata, obj.N5Info.Datatype);
        end
         
        function data = getIORegion(obj, regionStart, regionEnd)
            % In python indices
            regionStart = uint32(regionStart - 1);
            regionEnd = uint32(regionEnd);
            pydata = py.zarrAPI.getZarrRegion(obj.N5Obj, regionStart, regionEnd);
             
            data = cast(pydata, obj.N5Info.Datatype);
        end
        
    end
     
    % Write methods
    methods
        function openToWrite(obj, loc, info, level)
            assert(level==1)
            % TODO - map info.Datatype to corresponding N5 type.
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
            
            if ispc 
                if  numel(char(loc)) > 200
                    [pth, loc] = fileparts(loc);
                    cd(pth);
                end
            else
                if strcmp(loc{1}(1), '~')
                    homedir = getenv('HOME');
                    loc = sprintf("%s%s", homedir, loc{1}(2 : end));
                end
            end
            
            store = py.zarr.N5Store(loc);
            obj.N5Obj = py.zarr.open( ...
                pyargs('mode', 'w', ...
                'shape', uint32(info.Size), ...
                'chunks', uint32(info.IOBlockSize), ...
                'store', store, ...
                'compressor', py.numcodecs.Blosc(cname='zstd', clevel=1), ...
                'order', 'C', ...
                'dtype', dtype));
            obj.N5Info =  info;
        end
        
%         function openInParallelToAppend(obj, loc)
%             % assert(level==1)
%             % TODO - map info.Datatype to corresponding N5 type.
% %             obj.N5Obj = py.N5.open(loc,...
% %                 pyargs('mode', 'w',...
% %                 'shape', uint32(info.Size),...
% %                 'chunks', uint32(info.IOBlockSize),...
% %                 'dtype', 'u2')); % uint16 hardcoded
%             obj.N5Obj = py.N5.open(loc);
%             % obj.N5Info =  info;
%             obj.N5Info = obj.getInfo();
% 
%         end
        function openInParallelToAppend(obj, loc)
            obj.N5Obj = py.zarr.open(loc);
            obj.N5Info = obj.getInfo();
        end

        function setIOBlock(obj, ioBlockSub, level, data)
            assert(level==1)
            % In python indices
            if numel(ioBlockSub) == 2
                ioBlockSub = [ioBlockSub, 1];
            end
            regionStart = uint32((ioBlockSub-1).*obj.N5Info.IOBlockSize);
            regionEnd = uint32(double(regionStart) + obj.N5Info.IOBlockSize);
            % TODO - handle singleton dimensions, not sure how to convert
            % say a 4x1 data as 2D in python. Error looks like:
            %    Python Error: ValueError: parameter 'value': expected array with shape (4,1), got (4,)
            pydata = py.numpy.array(data);
            py.zarrAPI.setZarrRegion(obj.N5Obj, regionStart, regionEnd, pydata);
        end
        
        function setRegion(obj, regionStart, regionEnd, data)
            % In python indices
            regionStart = uint32(regionStart-1);
            regionEnd = uint32(regionEnd);
            % TODO - handle singleton dimensions, not sure how to convert
            % say a 4x1 data as 2D in python. Error looks like:
            %    Python Error: ValueError: parameter 'value': expected array with shape (4,1), got (4,)
            pydata = py.numpy.array(data);
            py.zarrAPI.setZarrRegion(obj.N5Obj, regionStart, regionEnd, pydata);
        end
        
        function setData(obj, data)
            % write the whole volume with multi-processing
            % In python indices
            % TODO - handle singleton dimensions, not sure how to convert
            % say a 4x1 data as 2D in python. Error looks like:
            %    Python Error: ValueError: parameter 'value': expected array with shape (4,1), got (4,)
            pydata = py.numpy.array(data);
            py.zarrAPI.setZarrData(obj.N5Obj, pydata);            
        end
    end
end
