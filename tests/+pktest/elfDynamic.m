function [symbolic, imports] = elfDynamic(fn)
% dynamic linking information of a 64-bit little-endian ELF file (Linux
% mex or shared library), read from its section headers, so it works on
% any OS.
%
% symbolic : true|false. It binds its own symbols first, before those of
%            libraries loaded earlier (linked with -Bsymbolic: DT_SYMBOLIC,
%            or DF_SYMBOLIC in DT_FLAGS).
%  imports : cell array. Names of the symbols it takes from other
%            libraries (undefined in its dynamic symbol table).

SHT_DYNAMIC = 6;
SHT_DYNSYM = 11;
DT_SYMBOLIC = 16;
DT_FLAGS = 30;
DF_SYMBOLIC = 2;

fid = fopen(fn, 'r', 'ieee-le');
if fid < 0
    error('Cannot open %s', fn);
end
cleanupObj = onCleanup(@() fclose(fid));

ident = fread(fid, [1, 16], 'uint8=>uint8');
if numel(ident) < 16 || ~isequal(ident(1 : 4), uint8([127, 'ELF'])) || ident(5) ~= 2 || ident(6) ~= 1
    error('%s is not a 64-bit little-endian ELF file', fn);
end
% e_shoff; e_shentsize and e_shnum
fseek(fid, 40, 'bof');
shoff = fread(fid, 1, 'uint64');
fseek(fid, 58, 'bof');
shentsize = fread(fid, 1, 'uint16');
shnum = fread(fid, 1, 'uint16');

% type, offset, size and link (the section of its strings) of each section
sections = zeros(shnum, 4);
for i = 1 : shnum
    fseek(fid, shoff + (i - 1) * shentsize + 4, 'bof');
    sections(i, 1) = fread(fid, 1, 'uint32');
    fseek(fid, shoff + (i - 1) * shentsize + 24, 'bof');
    sections(i, 2 : 3) = fread(fid, 2, 'uint64');
    sections(i, 4) = fread(fid, 1, 'uint32');
end

symbolic = false;
imports = {};
for i = find(sections(:, 1) == SHT_DYNAMIC)'
    % entries of tag and value
    d = reshape(typecast(readSection(i), 'uint64'), 2, []);
    symbolic = symbolic || any(d(1, :) == DT_SYMBOLIC) || any(bitand(d(2, d(1, :) == DT_FLAGS), DF_SYMBOLIC));
end
for i = find(sections(:, 1) == SHT_DYNSYM)'
    % 24-byte entries: name offset (uint32), info, other, section index
    % (uint16), value and size; undefined symbols have section index 0
    e = reshape(readSection(i), 24, []);
    nameOffsets = double(typecast(reshape(e(1 : 4, :), 1, []), 'uint32'));
    shndx = typecast(reshape(e(7 : 8, :), 1, []), 'uint16');
    strs = readSection(sections(i, 4) + 1);
    nuls = find(strs == 0);
    for o = nameOffsets(shndx == 0 & nameOffsets > 0)
        imports{end + 1} = char(strs(o + 1 : nuls(find(nuls > o, 1)) - 1));
    end
end

    function [data] = readSection(k)
        fseek(fid, sections(k, 2), 'bof');
        data = fread(fid, [1, sections(k, 3)], 'uint8=>uint8');
    end
end
