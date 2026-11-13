function [txt] = variant_mex_source(txt, mode, name)
% for the CI job exit-probe (build_variant_mex): the committed mex source txt
% changed for mode ('lock', 'atexit' or 'both'; see build_variant_mex). The
% mexLock block, unique in each source, is "#ifdef __APPLE__" with the comment
% "keep the mex file loaded"; the mex entry is one "void mexFunction(" at the
% start of a line. name is only for the error.

appleBlock = '#ifdef __APPLE__\r?\n\s*// keep the mex file loaded';
if numel(regexp(txt, appleBlock)) ~= 1 || numel(regexp(txt, '(?m)^void\s+mexFunction\s*\(')) ~= 1
    error('build_variant_mex:source', '%s is not of the expected shape.', name);
end
if any(strcmp(mode, {'atexit', 'both'}))
    % the handler, at file scope just before the mex entry (Windows only)
    handler = [newline, ...
        '#if defined(_WIN32)', newline, ...
        '#include <omp.h>', newline, ...
        '// end this mex file''s OpenMP worker pool before MATLAB unloads the mex', newline, ...
        '// file, so its OpenMP worker threads do not end under the Windows loader lock', newline, ...
        'static void petakit5d_end_omp_pool(void) { omp_pause_resource_all(omp_pause_hard); }', newline, ...
        '#endif', newline, newline];
    txt = regexprep(txt, '(?m)^(void\s+mexFunction\s*\()', [handler, '$1'], 'once');
    % register it, after the mexLock block (before the lock change below, so
    % the "#ifdef __APPLE__" anchor still matches)
    reg = [newline, ...
        '#if defined(_WIN32)', newline, ...
        '    mexAtExit(petakit5d_end_omp_pool);', newline, ...
        '#endif'];
    txt = regexprep(txt, '(?s)(#ifdef __APPLE__\r?\n\s*// keep the mex file loaded.*?#endif)', ['$1', reg], 'once');
end
if any(strcmp(mode, {'lock', 'both'}))
    % replace the opener only, keeping the comment ($1), as build_lockwin_mex did
    txt = regexprep(txt, '#ifdef __APPLE__(\r?\n\s*// keep the mex file loaded)', ...
        '#if defined(__APPLE__) || defined(_WIN32)$1', 'once');
end

end
