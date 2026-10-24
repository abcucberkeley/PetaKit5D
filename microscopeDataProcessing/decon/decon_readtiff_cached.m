function [im] = decon_readtiff_cached(fn)
% readtiff of a small file that decon reads again and again in a session (the
% generated PSF and the OMW back projector, for each batch of a large file or
% each file of a job): the last four files read are kept with their
% modification time and size, and a file is read again once they change.


persistent fns infos ims

if isempty(fns)
    fns = {};
    infos = {};
    ims = {};
end

dir_info = dir(fn);
if ~isscalar(dir_info) || dir_info.isdir
    im = readtiff(fn);
    return;
end
info = [dir_info.datenum, dir_info.bytes];
ind = find(strcmp(fns, fn), 1);
if ~isempty(ind) && isequal(infos{ind}, info)
    im = ims{ind};
    return;
end

im = readtiff(fn);
if ~isempty(ind)
    fns(ind) = [];
    infos(ind) = [];
    ims(ind) = [];
end
fns = [{fn}, fns(1 : min(end, 3))];
infos = [{info}, infos(1 : min(end, 3))];
ims = [{im}, ims(1 : min(end, 3))];

end
