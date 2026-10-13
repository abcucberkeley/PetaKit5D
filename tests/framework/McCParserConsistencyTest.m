classdef McCParserConsistencyTest < pktest.TestCase
    % the mcc parsers (mcc/parsers/*_parser.m) and mccMaster.m are generated
    % by mcc_autogen from mcc/mcc_function_list.txt and the input parsers of
    % the listed functions. If a function's parameters change and the
    % parsers are not regenerated, mcc mode (MATLAB Runtime) silently uses
    % stale defaults or rejects new parameters. These tests regenerate the
    % parsers in a temporary folder and compare them with the committed
    % files, without compiling anything.

    properties (Constant)
        % listed in mcc_function_list.txt, but their source is not in this
        % repository, so their parsers can't be regenerated
        MissingSourceFunctions = {'image_intensity_correction_block', 'replicateDataBlock'}
    end

    methods (Test, TestTags = {'Unit'})
        function functionListMatchesMccMasterCases(testCase)
            funcNames = testCase.listedFunctions();
            testCase.verifyEqual(numel(unique(funcNames)), numel(funcNames), 'Duplicate functions in the list.');
            cases = regexp(fileread(testCase.mccDir('parsers/mccMaster.m')), '^\s*case ''(\w+)''', 'tokens', 'lineanchors');
            cases = cellfun(@(x) x{1}, cases, 'unif', 0);
            testCase.verifyEqual(sort(cases(:)), sort(funcNames(:)));
        end

        function everyListedFunctionHasAParser(testCase)
            funcNames = testCase.listedFunctions();
            missing = funcNames(~cellfun(@(x) exist(testCase.mccDir(['parsers/', x, '_parser.m']), 'file') == 2, funcNames));
            testCase.verifyEmpty(missing, sprintf('No parser for: %s', strjoin(missing, ', ')));
        end

        function parsersMatchFunctionInputParsers(testCase)
            outPath = testCase.createTempFolder();
            [funcNames, dpntLists] = parse_mcc_function_list_file(testCase.mccDir('mcc_function_list.txt'));
            stale = {};
            for i = 1 : numel(funcNames)
                if any(strcmp(funcNames{i}, testCase.MissingSourceFunctions))
                    continue;
                end
                evalc('generate_function_mcc_parser_file(funcNames{i}, dpntLists{i}, outPath)');
                generated = testCase.readLines(sprintf('%s/%s_parser.m', outPath, funcNames{i}));
                committed = testCase.readLines(testCase.mccDir(['parsers/', funcNames{i}, '_parser.m']));
                if ~isequal(generated, committed)
                    n = min(numel(generated), numel(committed));
                    k = find(~strcmp(generated(1 : n), committed(1 : n)), 1);
                    if isempty(k)
                        k = n + 1;
                        stale{end + 1} = sprintf('%s: %d vs %d lines', funcNames{i}, numel(generated), numel(committed)); %#ok<AGROW>
                    else
                        stale{end + 1} = sprintf('%s line %d:\n    generated: %s\n    committed: %s', ...
                            funcNames{i}, k, generated{k}, committed{k}); %#ok<AGROW>
                    end
                end
            end
            testCase.verifyEmpty(stale, sprintf(['Parsers out of date with their functions (run mcc_autogen and ', ...
                'recompile):\n%s'], strjoin(stale, newline)));
        end

        function mccMasterMatchesRegenerated(testCase)
            outPath = testCase.createTempFolder();
            [funcNames, dpntLists] = parse_mcc_function_list_file(testCase.mccDir('mcc_function_list.txt'));
            committed = testCase.readLines(testCase.mccDir('parsers/mccMaster.m'));
            posArgNums = zeros(numel(funcNames), 1);
            paramNums = zeros(numel(funcNames), 1);
            for i = 1 : numel(funcNames)
                if any(strcmp(funcNames{i}, testCase.MissingSourceFunctions))
                    % no source: take the argument counts from the committed dispatch line
                    k = find(contains(committed, sprintf('''%s''', funcNames{i}))) + 1;
                    posArgNums(i) = numel(regexp(committed{k}, 'varargin\{\d+\}'));
                    paramNums(i) = contains(committed{k}, ':end}');
                    continue;
                end
                [~, ~, argTypes] = evalc_generate(funcNames{i}, dpntLists{i}, outPath);
                posArgNums(i) = sum(~strcmp(argTypes, 'parameter'));
                paramNums(i) = sum(strcmp(argTypes, 'parameter'));
            end
            evalc('generate_mcc_master_file(funcNames, posArgNums, paramNums, outPath)');
            testCase.verifyEqual(testCase.readLines([outPath, 'mccMaster.m']), committed, ...
                'mccMaster.m is out of date with the function list or the functions'' arguments.');
        end
    end

    methods (Test, TestTags = {'Unit', 'KnownIssue'})
        function everyListedFunctionHasSource(testCase)
            % KNOWN ISSUE: image_intensity_correction_block and
            % replicateDataBlock are listed in mcc/mcc_function_list.txt and
            % dispatched by mccMaster.m, but their source is not in this
            % repository, so compile_mccMaster can't include them and
            % mcc_autogen fails on them. Fix: add the sources or remove them
            % from the list and regenerate.
            funcNames = testCase.listedFunctions();
            missing = funcNames(cellfun(@(x) isempty(which(x)), funcNames));
            testCase.verifyEmpty(missing, sprintf('Listed functions without source: %s', strjoin(missing, ', ')));
        end

        function noOrphanParsers(testCase)
            % KNOWN ISSUE: mcc/parsers/multiple_region_spectra_encoding_block_parser.m
            % is not in the function list or mccMaster.m and its function
            % does not exist. Fix: delete the stale parser.
            dir_info = dir(testCase.mccDir('parsers/*_parser.m'));
            parsed = regexprep({dir_info.name}, '_parser\.m$', '');
            orphans = setdiff(parsed, testCase.listedFunctions());
            testCase.verifyEmpty(orphans, sprintf('Parsers for unlisted functions: %s', strjoin(orphans, ', ')));
        end
    end

    methods (Access = private)
        function [fn] = mccDir(~, name)
            fn = [strrep(pktest.repoRoot(), '\', '/'), '/mcc/', name];
        end

        function [funcNames] = listedFunctions(testCase)
            funcNames = parse_mcc_function_list_file(testCase.mccDir('mcc_function_list.txt'));
        end

        function [lines] = readLines(~, fn)
            % lines without carriage returns (CRLF checkouts on Windows)
            lines = readTextFile(fn);
            lines = strrep(lines(:), char(13), '');
        end
    end
end


function [funcLines, argLines, argTypes] = evalc_generate(funcName, dpntFuncNames, outPath)
% generate_function_mcc_parser_file without its console output

[~, funcLines, argLines, argTypes] = evalc('generate_function_mcc_parser_file(funcName, dpntFuncNames, outPath)');

end
