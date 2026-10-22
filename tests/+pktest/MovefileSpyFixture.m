classdef MovefileSpyFixture < matlab.unittest.fixtures.Fixture
    % wrap movefile in a spy that records each move and the files next to
    % its destination at that moment, then moves as usual. Tasks publish
    % their result by moving it into place, and the job framework counts a
    % task as done (and may cancel its job) once the result exists, so
    % everything else a task writes must exist before that move. The spy is
    % first on the path while the fixture is set up.

    properties (SetAccess = private)
        Folder
        CallsFile
    end

    properties (Access = private)
        WarningState
    end

    methods
        function setup(fixture)
            fixture.Folder = strrep(tempname, '\', '/');
            mkdir(fixture.Folder);
            % the calls are kept outside the spy's folder, so that the folder
            % on the path doesn't change
            fixture.CallsFile = [strrep(tempname, '\', '/'), '.mat'];
            lines = {
                'function varargout = movefile(varargin)'
                '% test spy written by pktest.MovefileSpyFixture'
                'if nargin >= 2'
                sprintf('    spyFn = ''%s'';', fixture.CallsFile)
                '    calls = struct(''src'', {}, ''dst'', {}, ''files'', {});'
                '    if exist(spyFn, ''file'')'
                '        a = load(spyFn);'
                '        calls = a.calls;'
                '    end'
                '    dst = pktest.MovefileSpyFixture.normalize(varargin{2});'
                '    calls(end + 1) = struct(''src'', pktest.MovefileSpyFixture.normalize(varargin{1}), ''dst'', dst, ...'
                '        ''files'', {pktest.MovefileSpyFixture.listFiles(fileparts(dst))});'
                '    save(spyFn, ''calls'');'
                'end'
                '[varargout{1 : nargout}] = builtin(''movefile'', varargin{:});'
                'end'};
            fid = fopen([fixture.Folder, '/movefile.m'], 'w');
            fprintf(fid, '%s\n', lines{:});
            fclose(fid);
            % MATLAB warns that the spy shadows a built-in each time it looks
            % through the path again
            fixture.WarningState = warning('off', 'MATLAB:dispatcher:nameConflict');
            addpath(fixture.Folder, '-begin');
            fixture.SetupDescription = 'Wrapped movefile in a spy';
        end

        function teardown(fixture)
            rmpath(fixture.Folder);
            clear movefile
            warning(fixture.WarningState);
            if exist(fixture.Folder, 'dir')
                rmdir(fixture.Folder, 's');
            end
            if exist(fixture.CallsFile, 'file')
                delete(fixture.CallsFile);
            end
        end

        function [calls] = getCalls(fixture)
            % recorded moves, each with src, dst and files (the files and
            % folders next to dst before the move, see listFiles)
            calls = struct('src', {}, 'dst', {}, 'files', {});
            if exist(fixture.CallsFile, 'file')
                a = load(fixture.CallsFile);
                calls = a.calls;
            end
        end

        function [files] = filesWhenMoved(fixture, dst)
            % the files next to dst (relative to its folder) when something
            % was last moved to dst; empty if nothing was
            calls = fixture.getCalls();
            ind = find(strcmp({calls.dst}, pktest.MovefileSpyFixture.normalize(dst)), 1, 'last');
            files = {};
            if ~isempty(ind)
                files = calls(ind).files;
            end
        end

        function [tf] = wasMoved(fixture, dst)
            % whether something was moved to dst
            calls = fixture.getCalls();
            tf = any(strcmp({calls.dst}, pktest.MovefileSpyFixture.normalize(dst)));
        end
    end

    methods (Static)
        function [files] = listFiles(folder)
            % files and folders under folder, as sorted paths relative to it,
            % without the contents of zarr folders
            files = cell(0, 1);
            if ~isfolder(folder)
                return;
            end
            d = dir([folder, '/**']);
            d = d(~ismember({d.name}, {'.', '..'}));
            top = dir(folder);
            root = strrep(top(1).folder, '\', '/');
            paths = strcat(strrep({d.folder}, '\', '/'), '/', {d.name});
            paths = extractAfter(paths, numel(root) + 1);
            parents = cellfun(@fileparts, paths, 'unif', 0);
            files = sort(paths(~contains(parents, '.zarr')))';
        end

        function [p] = normalize(p)
            % forward slashes, no repeated (except a leading //) or trailing
            % separators: the code under test builds paths like 'dir//name'
            p = regexprep(strrep(char(p), '\', '/'), '(?<=.)/+', '/');
            p = strip(p, 'right', '/');
        end
    end

    methods (Access = protected)
        function [tf] = isCompatible(~, ~)
            tf = true;
        end
    end
end
