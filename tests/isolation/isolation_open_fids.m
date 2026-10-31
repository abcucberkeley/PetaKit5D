function [fids] = isolation_open_fids()
% identifiers of the files open with fopen (openedFiles where it exists:
% fopen('all') warns that it will be removed, which would change lastwarn)


if exist('openedFiles', 'builtin') || exist('openedFiles', 'file')
    fids = openedFiles();
else
    fids = fopen('all');
end
fids = fids(:)';

end
