%makeSteerableMovie(img, M, sigma, varargin) generates a MP4 movie of a full rotation of the steerable filter response to the input image

% Francois Aguet, 01/31/2012

function makeSteerableMovie(img, M, sigma, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('Angles', 36, @isscalar);
ip.addParamValue('BorderCondition', 'mirror', @(x) any(strcmpi(x, {'mirror', 'periodic', 'replicate', 'zeros'})));
ip.addParamValue('FileName', 'SteerableMovie', @ischar);
ip.parse(varargin{:});

[~,~,~,fb] = steerableDetector(img, M, sigma, ip.Results.Angles, ip.Results.BorderCondition);


[~,~] = mkdir('frames');
intRange = [min(fb(:)) max(fb(:))]; % dynamic range should be same for all images
 
nt = size(fb, 3);
fmt = ['%0' num2str(ceil(log10(nt))) 'd'];

fprintf('Generating frames:     ');
for i = 1:nt
    frame = uint8(scaleContrast(fb(:,:,i), intRange)); 
    imwrite(frame, ['frames' filesep 'frame_' num2str(i, fmt) '.tif']);
    fprintf('\b\b\b\b%3d%%', round(100*i/nt));
end
fprintf('\n');


fr = '15';
cmd = ['ffmpeg -y -r ' fr ' -i ' 'frames' filesep 'frame_' fmt '.tif'...
    ' -b 20000k ' ip.Results.FileName '.mp4'];
system(cmd);

fprintf('\n');

