function [deps, installName, rpaths] = machoDylibs(fn)
% dynamic library dependencies of a 64-bit Mach-O file (macOS mex or
% dylib), read from its load commands, so it works on any OS. For
% universal (fat) files the arm64 slice is used.
%
%        deps : cell array. Libraries it loads (LC_LOAD_DYLIB, LC_LOAD_WEAK_DYLIB, LC_REEXPORT_DYLIB).
% installName : char. Its own install name (LC_ID_DYLIB), '' for mex files.
%      rpaths : cell array. Its LC_RPATH entries.

LC_LOAD_DYLIB = uint32(hex2dec('c'));
LC_ID_DYLIB = uint32(hex2dec('d'));
LC_LOAD_WEAK_DYLIB = uint32(hex2dec('80000018'));
LC_RPATH = uint32(hex2dec('8000001c'));
LC_REEXPORT_DYLIB = uint32(hex2dec('8000001f'));

fid = fopen(fn, 'r', 'ieee-le');
if fid < 0
    error('Cannot open %s', fn);
end
cleanupObj = onCleanup(@() fclose(fid));

% universal binary: big-endian header with one entry per architecture
start = 0;
if fread(fid, 1, 'uint32=>uint32', 0, 'ieee-be') == uint32(hex2dec('cafebabe'))
    nArch = fread(fid, 1, 'uint32', 0, 'ieee-be');
    for i = 1 : nArch
        % cputype, cpusubtype, offset, size, align
        arch = fread(fid, 5, 'uint32', 0, 'ieee-be');
        if arch(1) == hex2dec('0100000c')
            start = arch(3);
        end
    end
    if start == 0
        error('%s has no arm64 slice', fn);
    end
end

fseek(fid, start, 'bof');
magic = fread(fid, 1, 'uint32=>uint32');
if magic ~= uint32(hex2dec('feedfacf'))
    error('%s is not a 64-bit Mach-O file', fn);
end
% cputype, cpusubtype, filetype, ncmds, sizeofcmds, flags, reserved
header = fread(fid, 7, 'uint32=>uint32');
ncmds = double(header(4));

deps = {};
installName = '';
rpaths = {};
offset = start + 32;
for i = 1 : ncmds
    fseek(fid, offset, 'bof');
    cmd = fread(fid, 1, 'uint32=>uint32');
    cmdsize = double(fread(fid, 1, 'uint32=>uint32'));
    if any(cmd == [LC_LOAD_DYLIB, LC_ID_DYLIB, LC_LOAD_WEAK_DYLIB, LC_RPATH, LC_REEXPORT_DYLIB])
        strOffset = double(fread(fid, 1, 'uint32=>uint32'));
        fseek(fid, offset + strOffset, 'bof');
        % strings are relative to the load command
        str = fread(fid, cmdsize - strOffset, 'uint8=>char')';
        str = str(1 : find([str, char(0)] == char(0), 1) - 1);
        if cmd == LC_ID_DYLIB
            installName = str;
        elseif cmd == LC_RPATH
            rpaths{end + 1} = str;
        else
            deps{end + 1} = str;
        end
    end
    offset = offset + cmdsize;
end

end
