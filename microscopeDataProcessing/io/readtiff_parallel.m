%[s] = readtiff_parallel(filepath) loads a tiff file or stack using libtiff
% This function is ~1.5-2x faster than imread, useful for large stacks,
% and supports a wider range of TIFF formats (see below)
%
% Inputs:
%     filepath : path to TIFF file to read from
%
% 


function s = readtiff_parallel(filepath, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('filepath');
ip.addOptional('range', []);
ip.addOptional('info', [], @isstruct);
ip.addParameter('ShowWaitbar', false, @islogical);
ip.parse(filepath, varargin{:});
info = ip.Results.info;
range = ip.Results.range;

spmd
    warning('off', 'all'); % ignore unknown TIFF tags
end

if isempty(info)
    info = imfinfo(filepath);
end
if isempty(range)
    N = numel(info);
    range = 1:N;
else
    N = numel(range);
end

nx = info(1).Width;
ny = info(1).Height;


% Determine format and allocate arrays
isRGB = false;
if (strcmpi(info(1).ColorType, 'grayscale')||strcmpi(info(1).ColorType, 'indexed'))
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
        % s = zeros(ny,nx,N,cname);
        s = cell(N, 1);
    else % assume uint
        % s = zeros(ny,nx,N,['uint' num2str(info(1).BitDepth)]);
        s = cell(N, 1);
    end
elseif strcmpi(info(1).ColorType, 'truecolor') && info(1).BitDepth==24 && ...
        all(info(1).BitsPerSample==[8 8 8])
    s = cell(N,1);
    isRGB = true;
else
    error('Unsupported format');
end


% mf = memoize(@Tiff);

% t = Tiff(filepath, 'r');
if ~ip.Results.ShowWaitbar
    parfor i = 1:numel(range)
        t = Tiff(filepath, 'r');
    
        t.setDirectory(range(i));
        s{i} = t.read();
    end
%     p = gcp();
%     nworker = p.NumWorkers;
%     batchsize = ceil(numel(range) / nworker);
%     s = cell(nworker, 1);
%     parfor i = 1:nworker
%         t = Tiff(filepath, 'r');
%         range_i = range((i - 1) * batchsize + 1 : min(numel(range), i * batchsize));
%         for j = 1 : numel(range_i)
%             t.setDirectory(range_i(j));
%             s{i}{j} = t.read();
%         end
%     end
%     s = cat(2, s{:});
else
    [~,fname] = fileparts(filepath);
    h = waitbar(0, ['Loading ' fname]);
    for i = 1:numel(range)
        t.setDirectory(range(i));
        % s(:,:,i) = t.read();
        s{i} = t.read();
        waitbar(i/numel(range))
    end
    close(h);
end
% t.close();
s = cat(3, s{:});

% warning(w);
