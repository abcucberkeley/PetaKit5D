function fset = loadFigureSettings(mode)

if nargin<1
    mode = '';
end

% light blue
fset.cfB = hsv2rgb([5/9 0.3 1]);
fset.ceB = hsv2rgb([5/9 1 1]);

% dark blue
fset.cfB2 = hsv2rgb([0.6 0.3 1]);
fset.ceB2 = hsv2rgb([0.6 1 1]);


fset.cfR = hsv2rgb([0 0.3 1]);
fset.ceR = hsv2rgb([0 1 1]);

fset.cfG = hsv2rgb([1/3 0.3 1]);
fset.ceG = hsv2rgb([1/3 1 1]);

fset.cfY = hsv2rgb([1/7 0.3 1]);
fset.ceY = hsv2rgb([1/7 1 1]);

fset.cf0 = hsv2rgb([0 0 0.7]);
fset.ce0 = hsv2rgb([0 0 0.3]);

fset.categories = {'valid', 'faulty', 'cut', 'persistent',...
                        'valid', 'faulty', 'cut', 'persistent'};

fset.fontName = {'FontName', 'Helvetica'};

switch mode
    case 'print'
        fset.units = 'centimeters';
        fset.fPos = [10 10 8 5.5]; % centimeters
        fset.axPos = [1.5 1.5 6 3.5];
        fset.ifont = [fset.fontName, 'FontSize', 6];
        fset.tfont = [fset.fontName, 'FontSize', 7];
        fset.sfont = [fset.fontName, 'FontSize', 8];
        fset.lfont = [fset.fontName, 'FontSize', 10];
        fset.TickLength = [0.015 0.025];
        fset.axOpts = ['Layer', 'top', 'TickDir', 'out', 'LineWidth', 1, fset.sfont, 'TickLength', fset.TickLength];
        fset.fOpts = {'Units', fset.units, 'Position', [10 10 sum(fset.axPos([1 3]))+0.5 sum(fset.axPos([2 4]))+0.5],...
            'PaperPositionMode', 'auto', 'Color', 'w', 'InvertHardcopy', 'off'};
        fset.axOpts = [fset.axOpts 'Units', fset.units, 'Position', fset.axPos];
    otherwise
        fset.units = 'normalized';
        fset.fPos = get(0, 'DefaultFigurePosition');
        fset.axPos = get(0, 'DefaultAxesPosition');
        fset.ifont = [fset.fontName, 'FontSize', 10];
        fset.tfont = [fset.fontName, 'FontSize', 12];
        fset.sfont = [fset.fontName, 'FontSize', 12];
        fset.lfont = [fset.fontName, 'FontSize', 16];
        fset.TickLength = [0.01 0.025];
        fset.axOpts = ['Layer', 'top', 'TickDir', 'out', 'LineWidth', 1, fset.sfont];
        fset.fOpts = {'Position', get(0, 'DefaultFigurePosition'), 'Color', 'w'};
end


ce = [0 1 0; 1 1 0; 1 0.5 0; 1 0 0; 0 1 1; 0 0.5 1; 0 0 1; 0.5 0 1];
% ce = [0 1 0; 1 0 0; 1 0.5 0; 1 1 0; 0 1 1; 0.5 0 1; 0 0 1; 0 0.5 1];
cf = rgb2hsv(ce);
cf(:,2) = 0.3;
cf = hsv2rgb(cf);
fset.ceTrackClasses = ce;
fset.cfTrackClasses = cf;