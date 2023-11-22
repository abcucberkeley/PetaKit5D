function  [fig] = visualize_OTF_and_mask_outline(abs_OTF, OTF_mask, visualize)
% visualize OTF Mask on OTF


if nargin < 3
    visualize = true;
end

% abs_OTF = abs(fftshift(OTF));
sz = size(abs_OTF, 1 : 3);
csz = ceil((sz + 1) / 2);

abs_OTF_xy = abs_OTF(:, :, csz(3));
abs_OTF_yz = squeeze(abs_OTF(:, csz(2), :));
abs_OTF_xz = squeeze(abs_OTF(csz(1), :, :));

OTF_mask_xy = squeeze(OTF_mask(:, :, csz(3)));
OTF_mask_yz = squeeze(OTF_mask(:, csz(2), :));
OTF_mask_xz = squeeze(OTF_mask(csz(1), :, :));

OTF_xy_supp = OTF_mask_xy  & ~imerode(OTF_mask_xy, strel('disk', 1));
[y, x] = find(OTF_xy_supp);
omw_act_xy_coords = [y, x];

OTF_yz_supp = OTF_mask_yz  & ~imerode(OTF_mask_yz, strel('disk', 1));
[y, x] = find(OTF_yz_supp);
omw_act_yz_coords = [y, x];

OTF_xz_supp = OTF_mask_xz  & ~imerode(OTF_mask_xz, strel('disk', 1));
[y, x] = find(OTF_xz_supp);
omw_act_xz_coords = [y, x];


fontName = 'Helvetica';
fontSize = 14;

ws = [sz(1), sz(1), sz(2)];
ws = ws ./ sum(ws) * 0.85;
hs = [sz(2), sz(3), sz(3)];
hs = hs ./ max(hs(:)) * 0.7;
aps = [0.04, 0.2, ws(1), hs(1);
       ws(1) + 0.075, 0.2, ws(2), hs(2);
       ws(1) + ws(2) + 0.11, 0.2, ws(3), hs(3);
       ];

if sz(3) == 1
    aps = [0.04, 0.2, 0.8, 0.7];
end
if visualize
    figure('visible', 'on');
else
    figure('visible', 'off');
end

if sz(3) ~= 1
    set(gcf, 'Renderer', 'painters', 'Position', [1 1 1920 1080]);
else
    set(gcf, 'Renderer', 'painters', 'Position', [1 1 800 800]);
end
set(gcf, 'color', 'w')

h_1 = axes('position', aps(1, :));

log_abs_OTF_xy_crop = log10(abs_OTF_xy ./ max(abs_OTF_xy(:)))';
log_abs_OTF_yz_crop = log10(abs_OTF_yz ./ max(abs_OTF_yz(:)))';
log_abs_OTF_xz_crop = log10(abs_OTF_xz ./ max(abs_OTF_xz(:)))';

imagesc(log_abs_OTF_xy_crop)
colormap('hot');
clim([-3, 0])

hold on, 
s = scatter(omw_act_xy_coords(:, 1), omw_act_xy_coords(:, 2), '.');
% s.MarkerEdgeColor = [0.4940 0.1840 0.5560];
% s.MarkerFaceColor = [0.4940 0.1840 0.5560];

axis equal
xlim([1, size(log_abs_OTF_xy_crop, 2)])
ylim([1, size(log_abs_OTF_xy_crop, 1)])
set(gca, 'YDir', 'normal')
% xticks([1, 101, 201, 301, 401])
% xticklabels({'-2', '-1', '0', '1', '1'})
% yticks([1, 101, 201, 301, 401])
% yticklabels({'-1', '-0.5', '0', '0.5', '1'})

h_1.XAxis.FontSize = fontSize;
h_1.XAxis.FontName = fontName;
h_1.YAxis.FontSize = fontSize;
h_1.YAxis.FontName = fontName;
% xlabel('y (\mum)', 'FontName', fontName, 'FontSize', fontSize)
% xlabel('k_x/(4\pin/\lambda_{exc})', 'FontName', fontName, 'FontSize', fontSize)
% ylabel('k_y/(4\pin/\lambda_{exc})', 'FontName', fontName, 'FontSize', fontSize)
h_1.XTick = [];
h_1.YTick = [];
xlabel('k_x', 'FontName', fontName, 'FontSize', fontSize)
ylabel('k_y', 'FontName', fontName, 'FontSize', fontSize)

set(h_1, 'linewidth', 1);
h_1.TickDir = 'out';
h_1.TickLength = [0.0100, 0.0250];

if sz(3) == 1
    hc = colorbar;
    hc.Label.String = 'log_{10}(OTF Intensity)';
    hc.FontSize = fontSize;
    hc.FontName = fontName;
    
    hc.LineWidth = 1;
    hc.TickDirection = 'out';
    hc.TickLength = 0.01;
    
    set(gcf, 'InvertHardcopy', 'off');
    
    fig = gcf;
    return;
end


h_2 = axes('position', aps(2, :));
imagesc(log_abs_OTF_yz_crop)
colormap('hot');
clim([-3, 0])

hold on, 
s = scatter(omw_act_yz_coords(:, 1), omw_act_yz_coords(:, 2), '.');
% s.MarkerEdgeColor = [0.4940 0.1840 0.5560];
% s.MarkerFaceColor = [0.4940 0.1840 0.5560];

axis equal
xlim([1, size(log_abs_OTF_yz_crop, 2)])
ylim([1, size(log_abs_OTF_yz_crop, 1)])
set(gca, 'YDir', 'normal')
% xticks([1, 101, 201, 301, 401])
% xticklabels({'-1', '-0.5', '0', '0.5', '1'})
% yticks([1, 101, 201, 301, 401])
% yticklabels({'-1', '-0.5', '0', '0.5', '1'})
h_2.XTick = [];
h_2.YTick = [];
h_2.XAxis.FontSize = fontSize;
h_2.XAxis.FontName = fontName;
h_2.YAxis.FontSize = fontSize;
h_2.YAxis.FontName = fontName;
% xlabel('y (\mum)', 'FontName', fontName, 'FontSize', fontSize)
xlabel('k_x', 'FontName', fontName, 'FontSize', fontSize)
ylabel('k_z', 'FontName', fontName, 'FontSize', fontSize)

set(h_2, 'linewidth', 1);
h_2.TickDir = 'out';
h_2.TickLength = [0.0100, 0.0250];


h_3 = axes('position', aps(3, :));
imagesc(log_abs_OTF_xz_crop)
colormap('hot');
clim([-3, 0])

hold on, 
s = scatter(omw_act_xz_coords(:, 1), omw_act_xz_coords(:, 2), '.');
% s.MarkerEdgeColor = [0.4940 0.1840 0.5560];
% s.MarkerFaceColor = [0.4940 0.1840 0.5560];

axis equal
xlim([1, size(log_abs_OTF_xz_crop, 2)])
ylim([1, size(log_abs_OTF_xz_crop, 1)])
set(gca, 'YDir', 'normal')
% xticks([1, 101, 201, 301, 401])
% xticklabels({'-1', '-0.5', '0', '0.5', '1'})
% yticks([1, 101, 201, 301, 401])
% yticklabels({'-1', '-0.5', '0', '0.5', '1'})
h_3.XTick = [];
h_3.YTick = [];
h_3.XAxis.FontSize = fontSize;
h_3.XAxis.FontName = fontName;
h_3.YAxis.FontSize = fontSize;
h_3.YAxis.FontName = fontName;
% xlabel('y (\mum)', 'FontName', fontName, 'FontSize', fontSize)
xlabel('k_y', 'FontName', fontName, 'FontSize', fontSize)
ylabel('k_z', 'FontName', fontName, 'FontSize', fontSize)

set(h_3, 'linewidth', 1);
h_3.TickDir = 'out';
h_3.TickLength = [0.0100, 0.0250];

hc = colorbar;
hc.Label.String = 'log_{10}(OTF Intensity)';
hc.FontSize = fontSize;
hc.FontName = fontName;

hc.LineWidth = 1;
hc.TickDirection = 'out';
hc.TickLength = 0.01;

set(gcf, 'InvertHardcopy', 'off');

fig = gcf;

end
