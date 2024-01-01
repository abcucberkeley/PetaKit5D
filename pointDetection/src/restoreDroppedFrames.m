% Francois Aguet, 10/30/2013

function restoreDroppedFrames(data, frameID)

if numel(data)>1
    error('Input must be a single data set.');
end

if nargin<2
    frameID = '_stack';
end

nCh = numel(data.channels);
for c = 1:nCh
    
    % get frame #s and check for duplicates
    idx = str2double(regexpi(data.framePaths{c}, ['(?<=' frameID ')\d+'], 'match', 'once'));
    if min(idx)==0
        idx = idx+1;
    end
    rep = getMultiplicity(idx);
    cand = idx(rep>1);
    
    for i = 1:numel(cand)
        clist = data.framePaths{c}(idx==cand(i));
        if numel(clist)>2
            error('>1 frame to merge');
        end
        if numel(clist{1})<numel(clist{2})
            stackName = clist{1};
            frameName = clist{2};
        else
            stackName = clist{2};
            frameName = clist{1};
        end
        fprintf('Stack: %s\n', stackName);
        fprintf('Frame: %s\n', frameName);
       
        str = [];
        while ~any(strcmpi(str, {'y','n'}))
            str = input('Insert frame ? [y/n]: ', 's');
        end
        if strcmpi(str, 'y')
            stack = readtiff(stackName);
            frame = readtiff(frameName);
            
            % in some cases, >1 dropped but non-consecutive frames are
            % stored in a separate file
            ni = size(frame,3);
            for n = 1:ni
                stack = insertMissingFrame(stack, frame(:,:,n));
            end
            
            [~,fn,fe] = fileparts(stackName);
            stackName = [fn fe];
            [~,fn,fe] = fileparts(frameName);
            frameName = [fn fe];
            
            % move 
            bpath = [data.channels{c} 'errorFrames' filesep];
            [~,~] = mkdir(bpath);
            movefile([data.channels{c} frameName], [bpath frameName]);
            movefile([data.channels{c} stackName], [bpath stackName]);
            
            imwrite(stack(:,:,1), [data.channels{c} stackName], 'tif', 'compression' , 'none');
            for k = 2:size(stack,3)
                imwrite(stack(:,:,k), [data.channels{c} stackName], 'tif', 'compression' , 'none', 'writemode', 'append');
            end
        end
    end
end
