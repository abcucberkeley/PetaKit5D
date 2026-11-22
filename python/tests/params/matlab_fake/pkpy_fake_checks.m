function [out] = pkpy_fake_checks(x, varargin)
% A fake function for the probe helper's test (tests/params/test_parse_probe.py). Its checks
% return false or raise errors of their own: with an identifier, without one, one call further
% down, from a MathWorks function. It fails before its parse (an empty x) and after it (an
% error of its own, one of an inputParser of a function it calls, and one of its own parse in a
% call of itself: the same place, one frame further down).

first = x(1); %#ok<NASGU> (an empty x fails here, before the parse)
ip = inputParser;
ip.addRequired('x', @isnumeric);
ip.addParameter('Own', 1, @own_check);
ip.addParameter('Anon', 1, @(v) own_check(v));
ip.addParameter('Positive', 1, @(v) validateattributes(v, {'numeric'}, {'positive'}));
ip.addParameter('After', false, @islogical);
ip.addParameter('Nested', false, @islogical);
ip.addParameter('Recursive', false, @islogical);
ip.parse(x, varargin{:});
if ip.Results.After
    error('pkpyfake:afterParse', 'An error after the parse.');
end
if ip.Results.Nested
    nested_parse('text');
end
if ip.Results.Recursive
    pkpy_fake_checks('text');
end
out = 1;

end


function [tf] = own_check(v)
% false for a logical, its own errors for text, cells and structs, true otherwise
if ischar(v)
    error('pkpyfake:ownCheck', 'Text is refused.');
end
if iscell(v)
    error('A cell is refused.');
end
if isstruct(v)
    refuse_struct(v);
end
tf = ~islogical(v);
end


function refuse_struct(v)
error('pkpyfake:deeper', 'A struct of %d elements is refused.', numel(v));
end


function nested_parse(v)
ip = inputParser;
ip.addRequired('v', @isnumeric);
ip.parse(v);
end
