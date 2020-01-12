
function [m,isDone] = manualSegmentationTweakGUI(im,m,displayrange,isDone,ptsShow)
%MANUALSEGMENTATIONTWEAKGUI allows manual segmentation creation of masks or alteration of existing masksk
% [masks,isCompleted] = manualSegmentationTweakGUI(images,masks)
%
% This function allows the user to modify a set of masks manually to
% improve them. Based on Deepaks seed point selection GUI - thanks Deepak!!
%
%
%
% Instructions:
%
%   -Go to the frame you want to edit or create mask for (you can use mouse
%    scroll wheel to change frames).
%   -Select one of the options:
%
%               Add = add an area to the mask
%               Subtract = cut an area out of the mask
%               Restart = redraw this frame from scratch
%
%   -Select a drawing option:
%       Freehand = click and drag to select a region.
%       Polygon = click several times to create the vertices of a polygon.
%       Double click on first vertex to close polygon.
%
%
%   -Click GO or hit spacebar to start drawing your mask or mask correction
%
%   -Hit enter or click the "completed" box when you are done fixing a
%   frame
%
%   -When you are done with all the frames you want to fix, just close the
%   GUI
%
%   NOTE: To segment a cell which touches the image border, you must drag a
%   circle around it OUTSIDE the image area, or if using the polygon tool,
%   move a vertex outside of the image area.
%
% *****Keyboard shortcuts:*****
%
%   =For All the radio button options, just press the first letter
%   =Space - Go (start drawing on the mask)
%   =u - undo (only one step)
%   =m - toggle mask display
%   =enter - mark frame as completed
%   - OR + - Decrease/increase contrast by adjusting upper display limit
%   ( OR ) - Decrease/increase contrast by adjusting lower display limit
%
%
% Copyright (C) 2017, Danuser Lab - UTSouthwestern 
%
% This file is part of CMEAnalysis_Package.
% 
% CMEAnalysis_Package is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% CMEAnalysis_Package is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with CMEAnalysis_Package.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

%Hunter Elliott, 10/2012

% convert input to cell array of images if entered as 3D array
inputAsVol = false;
if ~iscell(im) 
    [ny,nx,nz] = size(im);
    im = squeeze(mat2cell(im, ny, nx, ones(nz,1)));
    inputAsVol = true;
end

if nargin < 5
    ptsShow = [];
end

if nargin < 4 || isempty(isDone)
    isDone = false(numel(im),1);
end

if nargin < 3 || isempty(displayrange)
    displayrange = cellfun(@(x) [min(x(:)) max(x(:))], im, 'unif', 0);
elseif ~iscell(displayrange)
    % duplicate range
    displayrange = mat2cell(repmat(displayrange, [numel(im) 1]), ones(numel(im),1), 2);
end

if nargin < 2 || isempty(m)
    m = cellfun(@(x) false(size(x)), im, 'unif', 0);
elseif ~iscell(m) % convert to cell array
    [ny,nx,nz] = size(m);
    m = squeeze(mat2cell(m, ny, nx, ones(nz,1)));
end

global guidata;

hMainFigure = fsFigure(.75);


% Create UI controls

% axis
guidata.ui.ah_img = axes( 'Position' , [ 0.001 , 0.2 , 0.7 , 0.7 ] , 'Visible' , 'off' );

% slice navigation controls
guidata.ui.pbh_dec = uicontrol(hMainFigure,'Style','pushbutton','String','<<',...
    'Units' , 'normalized' , 'Position',[0.20 0.1 0.05 0.05],...
    'Callback',{@pushFirstSlice_Callback});

guidata.ui.pbh_dec = uicontrol(hMainFigure,'Style','pushbutton','String','<',...
    'Units' , 'normalized' , 'Position',[0.25 0.1 0.05 0.05],...
    'Callback',{@pushdec_Callback});

guidata.ui.eth_sno = uicontrol(hMainFigure,'Style','edit',...
    'String','0',...
    'Units' , 'normalized' , 'Position',[0.30 0.1 0.1 0.05]);

guidata.ui.pbh_inc = uicontrol(hMainFigure,'Style','pushbutton','String','>',...
    'Units' , 'normalized' , 'Position',[0.40 0.1 0.05 0.05],...
    'Callback',{@pushinc_Callback});

guidata.ui.pbh_inc = uicontrol(hMainFigure,'Style','pushbutton','String','>>',...
    'Units' , 'normalized' , 'Position',[0.45 0.1 0.05 0.05],...
    'Callback',{@pushLastSlice_Callback});

% cursor point info controls
guidata.ui.eth_xloc = uicontrol(hMainFigure,'Style','edit',...
    'String','X: INV',...
    'Units' , 'normalized' , 'Position',[0.20 0.05 0.1 0.05]);

guidata.ui.eth_yloc = uicontrol(hMainFigure,'Style','edit',...
    'String','Y: INV',...
    'Units' , 'normalized' , 'Position',[0.30 0.05 0.1 0.05]);

guidata.ui.eth_Imval = uicontrol(hMainFigure,'Style','edit',...
    'String','I: INV',...
    'Units' , 'normalized' , 'Position',[0.40 0.05 0.1 0.05]);

% selection mode controls
guidata.ui.bgh_mode = uibuttongroup('visible','on', 'Units' , 'normalized' ,'Position',[0.71 0.2 0.2 0.2]);
guidata.ui_rbh_fgnd = uicontrol('Style','Radio','String','Add',...
    'Units' , 'normalized' ,'Position',[0.05 0.75 0.75 0.15],'parent',guidata.ui.bgh_mode,'HandleVisibility','off');
guidata.ui_rbh_bgnd = uicontrol('Style','Radio','String','Subtract',...
    'Units' , 'normalized' ,'Position',[0.05 0.50 0.75 0.15],'parent',guidata.ui.bgh_mode,'HandleVisibility','off');
guidata.ui_rbh_none = uicontrol('Style','Radio','String','Restart',...
    'Units' , 'normalized' ,'Position',[0.05 0.25 0.75 0.15],'parent',guidata.ui.bgh_mode,'HandleVisibility','off');

set( guidata.ui.bgh_mode , 'SelectedObject' , guidata.ui_rbh_fgnd );

% selection type controls
guidata.ui.sel_mode = uibuttongroup('visible','on', 'Units' , 'normalized' ,'Position',[0.71 0.5 0.2 0.2]);
guidata.ui_rbh2_fhan = uicontrol('Style','Radio','String','Freehand',...
    'Units' , 'normalized' ,'Position',[0.05 0.75 0.75 0.15],'parent',guidata.ui.sel_mode,'HandleVisibility','off');
guidata.ui_rbh2_poly = uicontrol('Style','Radio','String','Polygon',...
    'Units' , 'normalized' ,'Position',[0.05 0.50 0.75 0.15],'parent',guidata.ui.sel_mode,'HandleVisibility','off');

set( guidata.ui.sel_mode , 'SelectedObject' , guidata.ui_rbh2_fhan );


%Go button
guidata.ui_go = uicontrol('Style','pushbutton','String','Go',...
    'Units' , 'normalized' ,'Position',[0.71 0.8 0.2 0.1],'parent',hMainFigure,'Callback',{@pushGo_Callback});

%check box
guidata.ui_cb = uicontrol('Style','checkbox','String','Completed',...
    'Units' , 'normalized' ,'Position',[0.10 0.1 0.075 0.05],'parent',hMainFigure,'Callback',{@chkBox_Callback});

set(guidata.ui_cb,'Value',0)

%[0.20 0.1 0.05 0.05]
% set callbacks
set( hMainFigure , 'WindowScrollWheelFcn' , @FnSliceScroll_Callback );
%set( hMainFigure , 'WindowButtonDownFcn' , @FnMainFig_MouseButtonDownFunc );
%set( hMainFigure , 'WindowButtonMotionFcn' , @FnMainFig_MouseMotionFunc );
set( hMainFigure , 'WindowKeyPressFcn' , @FnKeyPress_Callback );


% guidata
guidata.im = im;
guidata.m = m;
guidata.prevm = m;
guidata.showMask = true;
guidata.isDone = isDone;
guidata.firstShow = true;%Stupid way to keep axis from resizing after initial
guidata.sliceno = 1;
guidata.displayrange = displayrange;
guidata.ptsShow = ptsShow;
guidata.ptsHan = [];
guidata.fgnd_seed_points = [];
guidata.bgnd_seed_points = [];

imsliceshow(guidata);
guidata.firstShow = false;


% wait until the window is closed
errCatch = 0;
try
    waitfor( hMainFigure );
catch
    errCatch = 1;
end

if errCatch == 0
    m = guidata.m;
    if inputAsVol % convert cell array back to 3D array
        m = cat(3, m{:});
    end
    isDone = guidata.isDone;
    clear guidata;
else
    clear guidata;
    error( 'Error: Unknown error occured while getting seed points from the user' );
end



function imsliceshow(guidata)

%Just show a blank mask so we can be lazy and still use deepaks display function
%and have consistant image display/contrast etc.
if guidata.showMask
    mShow = guidata.m{guidata.sliceno};
else
    mShow = false(size(guidata.m{guidata.sliceno}));
end

% FA: I commented the next 2 blocks such that the view area resets between
%     different images. This is required when input is a cell array of unequal sized
%     images

% xl = xlim(guidata.ui.ah_img);
% yl = ylim(guidata.ui.ah_img);

imHan = imshow(genImageMaskOverlay_loc(guidata.im{guidata.sliceno},...
    mShow,[0 1 0],.17,guidata.displayrange{guidata.sliceno}));

% if ~guidata.firstShow
%     xlim(guidata.ui.ah_img,xl)
%     ylim(guidata.ui.ah_img,yl)
% end

if ~isempty(guidata.ptsShow) && (isempty(guidata.ptsHan) || ~ishandle(guidata.ptsHan))
    hold on
    guidata.ptsHan = plot(guidata.ptsShow(:,1),guidata.ptsShow(:,2),'.r');
    hold off
end

set(guidata.ui.eth_sno,'String',sprintf('%d / %d' , guidata.sliceno , numel(guidata.im) ));
set(guidata.ui_cb,'Value',guidata.isDone(guidata.sliceno))

%Old method - show by transparency
%     imHan = imshow(guidata.im{guidata.sliceno},guidata.displayrange);
%     set(guidata.ui.eth_sno,'String',sprintf('%d / %d' , guidata.sliceno , size( guidata.im , 3 ) ));
%     set(imHan,'AlphaData',double(guidata.m{guidata.sliceno})+1);
%     set(imHan,'AlphaDataMapping','scaled')
%     alim(get(imHan,'Parent'),[0 2])
%     hold on;
%
%         if ~isempty( guidata.fgnd_seed_points )
%
%             fgnd_pt_ind = find( guidata.fgnd_seed_points( : , 3 ) == guidata.sliceno );
%             plot( guidata.fgnd_seed_points( fgnd_pt_ind , 1 ) , guidata.fgnd_seed_points( fgnd_pt_ind , 2 ) , 'g+' );
%
%         end
%
%         if ~isempty( guidata.bgnd_seed_points )
%
%             bgnd_pt_ind = find( guidata.bgnd_seed_points( : , 3 ) == guidata.sliceno );
%             plot( guidata.bgnd_seed_points( bgnd_pt_ind , 1 ) , guidata.bgnd_seed_points( bgnd_pt_ind , 2 ) , 'r+' );
%
%         end
%    hold off;



% First Slice
function pushFirstSlice_Callback(hSrc,eventguidata)
global guidata;

guidata.sliceno = 1;
imsliceshow(guidata);



% Last Slice
function pushLastSlice_Callback(hSrc,eventguidata)
global guidata;

guidata.sliceno = numel(guidata.im);
imsliceshow(guidata);



function pushdec_Callback(hSrc,eventguidata)
global guidata;

if(guidata.sliceno>1)
    guidata.sliceno = guidata.sliceno-1;
end
imsliceshow(guidata);



function pushinc_Callback(hSrc,eventguidata)

global guidata;

if(guidata.sliceno<numel(guidata.im))
    guidata.sliceno = guidata.sliceno+1;
end
imsliceshow(guidata);



function FnSliceScroll_Callback( hSrc , evnt )
global guidata;

if evnt.VerticalScrollCount > 0
    if(guidata.sliceno<numel(guidata.im))
        guidata.sliceno = guidata.sliceno+1;
    end    
elseif evnt.VerticalScrollCount < 0
    if(guidata.sliceno>1)
        guidata.sliceno = guidata.sliceno-1;
    end
end
imsliceshow(guidata);
UpdateCursorPointInfo(guidata);



% function FnMainFig_MouseButtonDownFunc( hSrc , evnt )
% 
% global guidata;
% 
% cp = get( gca , 'CurrentPoint' );
% 
% if IsPointInsideImage( cp(1,1:2) , guidata ) && strcmp( get(hSrc ,'SelectionType'),'normal' )
%     switch get( guidata.ui.bgh_mode , 'SelectedObject' )
%         
%         case guidata.ui_rbh_fgnd
%             guidata.fgnd_seed_points = [ guidata.fgnd_seed_points ; cp(1,1:2) guidata.sliceno ];
%             
%         case guidata.ui_rbh_bgnd
%             guidata.bgnd_seed_points = [ guidata.bgnd_seed_points ; cp(1,1:2) guidata.sliceno ];
%     end
% end
% imsliceshow(guidata);



% Update cursor point info -- xloc, yloc, int_val
function UpdateCursorPointInfo( guidata )

cp = get( gca , 'CurrentPoint' );

if IsPointInsideImage( cp(1,1:2) , guidata )
    set(guidata.ui.eth_xloc,'String' ,sprintf('X: %d / %d' , round( cp(1,1) ) , size(guidata.im{guidata.sliceno}, 2) ));
    set(guidata.ui.eth_yloc,'String' ,sprintf('Y: %d / %d' , round( cp(1,2) ) , size(guidata.im{guidata.sliceno}, 1) ));
    set(guidata.ui.eth_Imval,'String',sprintf('I: %.1f' , guidata.im{guidata.sliceno}( round(cp(1,2)), round(cp(1,1)) ) ));
else
    set(guidata.ui.eth_xloc,'String',sprintf('X: INV') );
    set(guidata.ui.eth_yloc,'String',sprintf('Y: INV') );
    set(guidata.ui.eth_Imval,'String',sprintf('I: INV') );
end



% function FnMainFig_MouseMotionFunc( hSrc , evnt )
% 
% global guidata;
% 
% cp = get( gca , 'CurrentPoint' );
% 
% if IsPointInsideImage( cp(1,1:2) , guidata )
%     set( hSrc ,'Pointer','crosshair');
% else
%     set( hSrc ,'Pointer','arrow');
% end
% imsliceshow(guidata);
% UpdateCursorPointInfo( guidata );



function [ blnInside ] = IsPointInsideImage( cp , guidata )
imsize = size(guidata.im{guidata.sliceno});
blnInside = all(cp<=imsize([2 1])) && all(cp>=[1 1]);



function pushGo_Callback(hSrc,eventguidata)

global guidata;

switch get( guidata.ui.sel_mode , 'SelectedObject' )
    case guidata.ui_rbh2_fhan
        fH = imfreehand(guidata.ui.ah_img);
    case guidata.ui_rbh2_poly
        fH = impoly(guidata.ui.ah_img);
end

if ~isempty(fH)
    currROI = fH.createMask;
    
    guidata.prevm = guidata.m;
    
    switch get( guidata.ui.bgh_mode , 'SelectedObject' )
        
        case guidata.ui_rbh_fgnd
            guidata.m{guidata.sliceno} = guidata.m{guidata.sliceno} | currROI;
            
        case guidata.ui_rbh_bgnd
            guidata.m{guidata.sliceno} = ...
                guidata.m{guidata.sliceno} ~= (currROI & guidata.m{guidata.sliceno});
            
        case guidata.ui_rbh_none            
            guidata.m{guidata.sliceno} = currROI;
    end
end
imsliceshow(guidata);



function FnKeyPress_Callback(hSrc,eventguidata)
global guidata;

switch eventguidata.Key
    case 'space'
        %Call the go button function
        pushGo_Callback(hSrc,eventguidata)
        
    case 'a'
        set( guidata.ui.bgh_mode , 'SelectedObject' , guidata.ui_rbh_fgnd);
        
    case 's'
        set( guidata.ui.bgh_mode , 'SelectedObject' , guidata.ui_rbh_bgnd);
        
    case 'r'
        set( guidata.ui.bgh_mode , 'SelectedObject' , guidata.ui_rbh_none);
        
    case 'f'
        set( guidata.ui.sel_mode , 'SelectedObject' , guidata.ui_rbh2_fhan);
        
    case 'p'
        set( guidata.ui.sel_mode , 'SelectedObject' , guidata.ui_rbh2_poly);
        
    case 'equal'
        guidata.displayrange{guidata.sliceno} = guidata.displayrange{guidata.sliceno} - [0 100];
        
    case 'hyphen'
        guidata.displayrange{guidata.sliceno} = guidata.displayrange{guidata.sliceno} + [0 100];
        
    case '0'
        if ~isempty(eventguidata.Character)
            guidata.displayrange{guidata.sliceno} = guidata.displayrange{guidata.sliceno} - [100 0];
        end
    case '9'
        guidata.displayrange{guidata.sliceno} = guidata.displayrange{guidata.sliceno} + [100 0];
        
    case 'm'
        guidata.showMask = ~guidata.showMask;
        
    case 'return'
        chkBox_Callback(hSrc,eventguidata)
        
    case 'u' % in case of regrets...
        tmp = guidata.m;
        guidata.m = guidata.prevm;
        guidata.prevm = tmp;
end
imsliceshow(guidata);



function chkBox_Callback(hSrc,eventguidata)
global guidata
guidata.isDone(guidata.sliceno) = ~guidata.isDone(guidata.sliceno);



function imRGB = genImageMaskOverlay_loc(im, mask, maskColor, maskAlpha, displayRange)

imr = im2uint8( mat2gray( im ,double(displayRange)) );
img = imr;
imb = imr;
mask = mask > 0;
maxVal = 255;

imr(mask) = uint8( double( (1 - maskAlpha) * imr(mask) ) + maxVal * maskAlpha * maskColor(1) );
img(mask) = uint8( double( (1 - maskAlpha) * img(mask) ) + maxVal * maskAlpha * maskColor(2) );
imb(mask) = uint8( double( (1 - maskAlpha) * imb(mask) ) + maxVal * maskAlpha * maskColor(3) );

imRGB = cat(3, imr, img, imb );
