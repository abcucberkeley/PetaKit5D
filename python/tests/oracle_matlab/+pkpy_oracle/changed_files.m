function [files, deleted] = changed_files(a, b)
% what a call wrote, from snapshots before (a) and after (b): new or changed
% files (a zarr folder as a whole) and new empty folders; and the files it
% removed. Same rules as petakit5d.testing.cases.changed_files.


files = {};
deleted = {};
parents = cellfun(@parent_of, b.rel, 'unif', 0);
[inA, ia] = ismember(b.rel, a.rel);
for i = 1 : numel(b.rel)
    rel = b.rel{i};
    z = zarr_root(rel);
    if b.isdir(i)
        if strcmp(z, rel)
            if ~inA(i)
                files{end + 1} = rel; %#ok<AGROW>
            end
        elseif isempty(z) && ~inA(i) && ~any(strcmp(parents, rel))
            files{end + 1} = rel; %#ok<AGROW>
        end
        continue;
    end
    if ~inA(i) || a.isdir(ia(i)) || a.bytes(ia(i)) ~= b.bytes(i) || a.datenum(ia(i)) ~= b.datenum(i)
        if isempty(z)
            files{end + 1} = rel; %#ok<AGROW>
        else
            files{end + 1} = z; %#ok<AGROW>
        end
    end
end
inB = ismember(a.rel, b.rel);
for i = reshape(find(~inB), 1, [])
    rel = a.rel{i};
    z = zarr_root(rel);
    if ~isempty(z)
        if any(strcmp(b.rel, z))
            files{end + 1} = z; %#ok<AGROW>
        elseif strcmp(z, rel)
            deleted{end + 1} = rel; %#ok<AGROW>
        end
    elseif ~a.isdir(i)
        deleted{end + 1} = rel; %#ok<AGROW>
    end
end
files = reshape(unique(files), 1, []);
deleted = reshape(unique(deleted), 1, []);

end


function [p] = parent_of(rel)

k = find(rel == '/', 1, 'last');
if isempty(k)
    p = '';
else
    p = rel(1 : k - 1);
end

end


function [z] = zarr_root(rel)

z = '';
parts = strsplit(rel, '/');
for i = 1 : numel(parts)
    if endsWith(lower(parts{i}), '.zarr')
        z = strjoin(parts(1 : i), '/');
        return;
    end
end

end
