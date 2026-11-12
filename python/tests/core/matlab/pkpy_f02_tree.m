function root = pkpy_f02_tree(spec)
% make the file tree of a spec (cell of 'f|path', 'd|path', 'l|path|target') in a new temp folder
root = tempname;
mkdir(root);
for i = 1:numel(spec)
    t = strsplit(spec{i}, '|');
    p = [root, '/', t{2}];
    switch t{1}
        case 'd'
            mkdir(p);
        case 'f'
            fid = fopen(p, 'w');
            fwrite(fid, uint8(t{2}));
            fclose(fid);
        case 'l'
            system(sprintf('ln -s ''%s'' ''%s''', t{3}, p));
    end
end
end
