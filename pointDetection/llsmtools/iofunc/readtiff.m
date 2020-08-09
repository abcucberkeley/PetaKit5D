%[s] = readtiff(filepath) loads a tiff file or stack using libtiff
% This function is ~1.5-2x faster than imread, useful for large stacks,
% and supports a wider range of TIFF formats (see below)
%
% Inputs:
%     filepath : path to TIFF file to read from
%
% Optional inputs
%        range : range of pages to read from multi-page TIFF (stack)
%         info : output from imfinfo for filepath
%
% Options ('specifier', value)
%   'ShowWaitbar' : true|{false} displays a progess bar while loading
%
%
% Supported formats: unsigned integer, signed integer, single, and double TIFFs
%                    as well as 8-bit/channel RGB stacks (returned as a cell array)

% Francois Aguet, 05/21/2013

function s = readtiff(filepath, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('filepath');
ip.addOptional('range', []);
ip.addOptional('info', [], @isstruct);
ip.addParamValue('ShowWaitbar', false, @islogical);
ip.parse(filepath, varargin{:});
info = ip.Results.info;
range = ip.Results.range;
if isempty(info)
    info = imfinfo(filepath);
end
if isempty(range)
    N = numel(info);
    range = 1:N;
else
    N = numel(range);
end

w = warning('off', 'all'); % ignore unknown TIFF tags

nx = info(1).Width;
ny = info(1).Height;


% Determine format and allocate arrays
isRGB = false;
if strcmpi(info(1).ColorType, 'grayscale')
    if isfield(info, 'SampleFormat')
        switch info(1).SampleFormat
            case 'Unsigned integer'
                cname = ['uint' num2str(info(1).BitDepth)];
            case 'IEEE floating point'
                if info(1).BitDepth==32
                    cname = 'single';
                else
                    cname = 'double';
                end
            case 'Two''s complement signed integer'
                cname = ['int' num2str(info(1).BitDepth)];
            otherwise
                error('Unsupported format.');
        end
        s = zeros(ny,nx,N,cname);
    else % assume uint
        s = zeros(ny,nx,N,['uint' num2str(info(1).BitDepth)]);
    end
elseif strcmpi(info(1).ColorType, 'truecolor') && info(1).BitDepth==24 && ...
        all(info(1).BitsPerSample==[8 8 8])
    s = cell(N,1);
    isRGB = true;
else
    error('Unsupported format');
end


t = Tiff(filepath, 'r');
if ~ip.Results.ShowWaitbar
    if ~isRGB
        for i = 1:numel(range)
            t.setDirectory(range(i));
            s(:,:,i) = t.read();
        end
    else
        for i = 1:numel(range)
            t.setDirectory(range(i));
            s{i} = t.read();
        end
    end
else
    [~,fname] = fileparts(filepath);
    h = waitbar(0, ['Loading ' fname]);
    for i = 1:numel(range)
        t.setDirectory(range(i));
        s(:,:,i) = t.read();
        waitbar(i/numel(range))
    end
    close(h);
end
t.close();

warning(w);
