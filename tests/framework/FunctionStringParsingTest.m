classdef FunctionStringParsingTest < pktest.TestCase
    % separate_function_string_variables and
    % convert_function_string_to_mcc_string, which turn the task strings of
    % the job framework into command line arguments for the compiled (mcc)
    % mccMaster in mcc mode

    properties (TestParameter)
        % task string and its expected arguments (raw text, quotes included)
        parsedCase = struct( ...
            'quotedNumberMatrixCell', {{'f(''a'',1,''b'',[1,2],{''c'',''d''})', {'''a''', '1', '''b''', '[1,2]', '{''c'',''d''}'}}}, ...
            ... % 73ef511 (error) and 5e38703 (duplicated first argument) regressions
            'unquotedFirstArgument', {{'f(5,''a'')', {'5', '''a'''}}}, ...
            'emptyStringArgument', {{'f(''a'','''',3)', {'''a''', '''''', '3'}}}, ...
            'emptyStringInCell', {{'f(''a'',{''''},''b'')', {'''a''', '{''''}', '''b'''}}}, ...
            'emptyStringLast', {{'f(''a'',''b'','''')', {'''a''', '''b''', ''''''}}}, ...
            'nestedQuotedFunctionString', {{'f(''x'',''fs'',''g(''''y'''',1)'')', {'''x''', '''fs''', '''g(''''y'''',1)'''}}}, ...
            'nameValuePairs', {{'XR_crop_frame(''/d/a.tif'',''/d/c/a.tif'',[1,1,1,8,8,8],''pad'',false,''padValue'',0)', ...
                {'''/d/a.tif''', '''/d/c/a.tif''', '[1,1,1,8,8,8]', '''pad''', 'false', '''padValue''', '0'}}})
    end

    methods (Test, TestTags = {'Unit'})
        function separatesFunctionNameAndArguments(testCase, parsedCase)
            [funcName, args] = separate_function_string_variables(parsedCase{1});
            testCase.verifyEqual(funcName, regexp(parsedCase{1}, '^\w+', 'match', 'once'));
            testCase.verifyEqual(args, parsedCase{2});
        end

        function unmixBlockTaskString(testCase)
            % shape of the task strings built by XR_unmix_channels_zarr; the
            % unquoted first argument was duplicated before 5e38703
            str = ['unmix_channels_block(3,{''/d/a.zarr'',''/d/b.zarr''},''/d/tmp.zarr'',[1,0.5],', ...
                '''/d/blocks_3_3.mat'',[1,1,1,8,8,8],''uuid'',''u1'',''debug'',false)'];
            [funcName, args] = separate_function_string_variables(str);
            testCase.verifyEqual(funcName, 'unmix_channels_block');
            testCase.verifyEqual(args, {'3', '{''/d/a.zarr'',''/d/b.zarr''}', '''/d/tmp.zarr''', '[1,0.5]', ...
                '''/d/blocks_3_3.mat''', '[1,1,1,8,8,8]', '''uuid''', '''u1''', '''debug''', 'false'});
        end

        function convertsToMccCommandLineArguments(testCase)
            [funcName, varStr] = convert_function_string_to_mcc_string('f(''a'',[1,2],{''c'',''d''},true,[1;2])');
            testCase.verifyEqual(funcName, 'f');
            % cells and matrices with ';' are quoted for the shell; on
            % Windows single quotes become double quotes
            if ispc
                expected = """a"" [1,2] ""{""c"",""d""}"" true ""[1;2]""";
            else
                expected = "'a' [1,2] ""{'c','d'}"" true ""[1;2]""";
            end
            testCase.verifyEqual(varStr, expected);
        end
    end

    methods (Test, TestTags = {'Unit', 'KnownIssue'})
        function splitsLeadingUnquotedArguments(testCase)
            % KNOWN ISSUE: separate_function_string_variables.m:97-101 takes
            % the whole text before the first quoted/bracketed argument as
            % one argument, so 'f(1,2,''a'')' gives {'1,2', '''a'''}. Fix:
            % strsplit that leading text on ',' as the between-pairs branch
            % does.
            [~, args] = separate_function_string_variables('f(1,2,''a'')');
            testCase.verifyEqual(args, {'1', '2', '''a'''});
        end

        function keepsArgumentsWithoutQuotesOrBrackets(testCase)
            % KNOWN ISSUE: with no quoted/bracketed argument at all
            % (p_pairs is empty), separate_function_string_variables.m:97
            % takes var_str(s : max(s, min([]) - 2)) and loses the
            % arguments: 'f(1,2)' gives a single empty argument. Fix:
            % handle isempty(p_pairs) by strsplit(var_str, ',').
            [~, args] = separate_function_string_variables('f(1,2)');
            testCase.verifyEqual(args, {'1', '2'});
        end

        function acceptsEmptyStringAsFirstArgument(testCase)
            % KNOWN ISSUE: '''''' (an empty char) is only recognized when it
            % follows a comma or is inside {} (patterns ,'', {''} and ,''$
            % at separate_function_string_variables.m:31-33), so an empty
            % first argument errors with "'' are not fully paired".
            [~, args] = separate_function_string_variables('f('''',''a'')');
            testCase.verifyEqual(args, {'''''', '''a'''});
        end
    end
end
