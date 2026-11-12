"""Reading MATLAB sources: comments, continuations, strings, inputParser and arguments blocks."""

from __future__ import annotations

from petakit5d.params.matlab_source import (
    clean_source,
    demo_literals,
    extract_arguments_block,
    extract_inputparser,
    extract_parser_conversions,
    mcc_mode,
    split_args,
)

FUNC = """\
function [] = XR_fake(dataPaths, varargin)
% Required inputs:
%           dataPaths : char or cell array. Folders.
% Parameters:
%           pad : true|false (default: false). Pad it.

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('pad', false, @islogical); % a comment with 'quotes'
ip.addParameter('names', {'a, b', 'c''d'}, ...
    @iscell);
%{
ip.addParameter('hidden', 1, @isnumeric);
%}
ip.addParameter('box', [1 2; 3 4]', @(x) isnumeric(x) && numel(x) == 4);
ip.parse(dataPaths, varargin{:});
end

function local_helper()
ip = inputParser;
ip.addParameter('notMine', 1);
end
"""

ARGS = """\
function data = readfake(filepath, options)
arguments
    filepath char
    options.range (1, :) {mustBeNumeric} = []
    options.mode char {mustBeMember(options.mode, {'a', 'b'})} = 'a'  % note
end
data = 1;
end
"""


def test_clean_source_keeps_strings():
    out = clean_source("x = 'a % not a comment'; % comment\ny = [1, ...\n 2]; z = a';\n")
    assert "'a % not a comment'" in out and "comment\n" not in out
    assert "..." not in out and "[1," in out and "2]" in out and "a'" in out
    assert out.count("\n") == 2  # the continuation joined two lines


def test_split_args():
    assert split_args("'a', {1, 2}, @(x) f(x, 2), 'it''s, ok'") == [
        "'a'",
        "{1, 2}",
        "@(x) f(x, 2)",
        "'it''s, ok'",
    ]


def test_extract_inputparser(tmp_path):
    f = tmp_path / "XR_fake.m"
    f.write_text(FUNC)
    fn = extract_inputparser(f)
    assert fn.name == "XR_fake" and fn.args == ["dataPaths", "varargin"]
    assert [(p.kind, p.name) for p in fn.params] == [
        ("required", "dataPaths"),
        ("parameter", "pad"),
        ("parameter", "names"),
        ("parameter", "box"),
    ]
    assert fn.params[2].default == "{'a, b', 'c''d'}"
    assert fn.params[2].validator == "@iscell"
    assert fn.params[3].default == "[1 2; 3 4]'"
    assert fn.params[3].validator == "@(x) isnumeric(x) && numel(x) == 4"
    assert fn.case_sensitive is False and fn.partial_matching is True


def test_extract_arguments_block(tmp_path):
    f = tmp_path / "readfake.m"
    f.write_text(ARGS)
    fn = extract_arguments_block(f)
    assert [(p.kind, p.name, p.default, p.validator) for p in fn.params] == [
        ("required", "filepath", None, "char"),
        ("parameter", "range", "[]", "(1, :) {mustBeNumeric}"),
        ("parameter", "mode", "'a'", "char {mustBeMember(options.mode, {'a', 'b'})}"),
    ]
    g = tmp_path / "plain.m"
    g.write_text("function out = plain(a, b, varargin)\nout = a;\nend\n")
    assert [(p.kind, p.name) for p in extract_arguments_block(g).params] == [
        ("positional", "a"),
        ("positional", "b"),
    ]


def test_parser_conversions(tmp_path):
    f = tmp_path / "XR_fake_parser.m"
    f.write_text(
        FUNC.replace("XR_fake", "XR_fake_parser")
        + "\nif ischar(pad)\n    pad = str2num(pad);\nend\n"
        "if ischar(names) && ~isempty(names) && strcmp(names(1), '{')\n"
        "    names = eval(names);\nend\n"
    )
    conv = extract_parser_conversions(f)["conversions"]
    assert conv == {"pad": "str2num", "names": "eval:{"}


def test_mcc_mode_rules():
    assert mcc_mode("@ischar") == "none"
    assert mcc_mode("@isnumeric") == "str2num"
    assert mcc_mode("@islogical") == "str2num"
    assert mcc_mode("@(x) ischar(x) || iscell(x)") == "eval:{"
    assert mcc_mode("@isempty") == "eval:[{"
    assert mcc_mode("@(x) isempty(x) || iscell(x)") == "eval:@[{"
    assert mcc_mode("@(x) isnumeric(x) || ischar(x)") == "none"
    assert mcc_mode("@(x) isempty(x) || isnumeric(x)") == "str2num"
    assert mcc_mode("@isstruct") == "eval:["
    assert mcc_mode("@(x) isempty(x) || isa(x, 'function_handle')") == "eval:@[{"
    assert mcc_mode(None) == "none"


def test_demo_literals(tmp_path):
    f = tmp_path / "demo_fake.m"
    f.write_text(
        "dz = 0.3; % step\nChannelPatterns = {'CamB_ch0', 'CamA_ch0'};\n"
        "XR_fake(dataPath, 'Reverse', true, 'dz', dz, Save16bit=[true, true]);\n"
    )
    items = [(i["kind"], i["where"], i["text"]) for i in demo_literals(f)]
    assert ("assign", "dz", "0.3") in items
    assert ("assign", "ChannelPatterns", "{'CamB_ch0', 'CamA_ch0'}") in items
    assert ("arg", "XR_fake", "true") in items
    assert ("arg", "XR_fake:Save16bit", "[true, true]") in items
