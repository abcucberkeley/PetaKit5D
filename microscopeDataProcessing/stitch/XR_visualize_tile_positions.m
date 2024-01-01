function [] = XR_visualize_tile_positions(imageListFileName, varargin)
% function to visualize tile positions in order to quickly inspect whether their
% positions are correct. 
% 
% 
% Author: Xiongtao Ruan


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('imageListFileName', @isstr);
% ip.addParameter('Overwrite', true, @islogical);

ip.parse(imageListFileName, varargin{:});


% read image list csv file and parse the information
T = readtable(imageListFileName, 'Delimiter','comma');
% t_column_name = t.Properties.VariableNames;
fn = T.Filename;

specifyCam = true;
if all(~cellfun(@isempty, regexp(fn, '_Cam\w_ch', 'match')))
    expression = 'Scan_Iter_(?<Iter>\d+)_Cam(?<Cam>\w+)_ch(?<ch>\d+)_CAM1_stack(?<stack>\d+)_(?<laser>\d+)nm_(?<abstime>\d+)msec_(?<fpgatime>\d+)msecAbs_(?<x>\d+)x_(?<y>\d+)y_(?<z>\d+)z_(?<t>\d+)t.tif';
elseif all(~cellfun(@isempty, regexp(fn, '_ch[0-9]_', 'match')))
    expression = 'Scan_Iter_(?<Iter>\d+)_ch(?<ch>\d+)_CAM1_stack(?<stack>\d+)_(?<laser>\d+)nm_(?<abstime>\d+)msec_(?<fpgatime>\d+)msecAbs_(?<x>\d+)x_(?<y>\d+)y_(?<z>\d+)z_(?<t>\d+)t.tif';
    specifyCam = false;
end
tmp = regexpi(fn, expression, 'names');

t = T;
for f = 1:numel(tmp)
    t.Iter(f) = str2double(tmp{f}.Iter);
    if specifyCam
        t.camera(f) = (tmp{f}.Cam);
    else
        % use N to represent cam if it is not contained in the filename
        t.camera(f) = 'N';
    end
    t.ch(f) = str2double(tmp{f}.ch);
    t.stack(f) = str2double(tmp{f}.stack);
    t.laser(f) = str2double(tmp{f}.laser);
    t.abstime(f) = str2double(tmp{f}.abstime);
    t.fpgatime(f) = str2double(tmp{f}.fpgatime);
    t.x(f) = str2double(tmp{f}.x);
    t.y(f) = str2double(tmp{f}.y);
    t.z(f) = str2double(tmp{f}.z);
    t.t(f) = str2double(tmp{f}.t);
end


% by default, visualize the first channel and first time point
t_1 = t(t.Iter == 0 & t.ch == t.ch(1) & t.camera == t.camera(1), :);

xyz = [t_1.StageX_um_, t_1.StageY_um_, t_1.StageZ_um_];

% make scatter plot
figure('PaperPositionMode', 'auto', 'Color', 'w', 'InvertHardcopy', 'off');
set(gcf, 'Position', [1, 1, 1200, 1000]);

axes('Position', [0.15, 0.15, 0.7, 0.7]);
scatter3(xyz(:, 1) + 7, xyz(:, 2) + 7, xyz(:, 3), 50, '*');
for i = 1 : size(xyz, 1)
    txt = sprintf('%03dx_%03dy_%03dz', t_1.x(i), t_1.y(i), t_1.z(i));
    text(xyz(i, 1), xyz(i, 2), xyz(i, 3), txt, 'Interpreter','none');
end

xlabel('x');
ylabel('y');
zlabel('z');
view([-45, 45]);
title('Tile positions');
set(findobj(gca, '-property', 'fontsize'), 'fontsize', 18)


figure_location = fileparts(imageListFileName);
figure_filename = sprintf('%s/visualize_tile_positions.eps', figure_location);
f0 = gcf();
print(f0, '-painters','-depsc', '-loose', figure_filename);


end

