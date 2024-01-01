function [res] = XR_FSC_resolution(fsc, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('fsc', @isstruct);
ip.addParameter('xyPixelSize', 0.108, @isnumeric);
ip.addParameter('dz', 0.108, @isnumeric);
ip.addParameter('resThreshMethod', 'fixed', @ischar);
ip.addParameter('resThresh', 0.2, @isnumeric);
ip.addParameter('debug', false, @islogical);

ip.parse(fsc, varargin{:});

pr = ip.Results;
xyPixelSize = pr.xyPixelSize;
dz = pr.dz;
resThreshMethod = pr.resThreshMethod;
resThresh = pr.resThresh;

frc_mat = fsc.fsc;
r = fsc.r;
thetas = fsc.thetas;
npoints = fsc.npoints;

% use threshold 0.2
switch resThreshMethod
    case 'fixed'
        thrsh = resThresh;
        points = zeros(numel(thetas), 1);
        for i = 1 : numel(thetas)
            frc_i = frc_mat(:, i);
            r_i = r;

            % fit with smoothing spline
            frc_s_i = fit(r_i, frc_i, 'smoothingSpline');

            rs_i = r_i(1) : 0.0001 : 1;
            ind = find(frc_s_i(rs_i) < thrsh, 1, 'first');
            if isempty(ind)
                rs_i = r_i(1) : 0.0001 : sqrt(2);
                ind = find(frc_s_i(rs_i) < thrsh, 1, 'first');
                if isempty(ind)
                    frc_s_i = fit(r_i, frc_i, 'poly2');
                    ind = find(frc_s_i(rs_i) < thrsh, 1, 'first');
                end
            end
            if isempty(ind)
                points(i) = sqrt(2);
            else 
                points(i) = rs_i(ind);
            end
        end
    case {'one-bit', 'half-bit'}
        points = zeros(numel(thetas), 1);
        for i = 1 : numel(thetas)
            frc_i = frc_mat(:, i);
            r_i = r;

            % fit with smoothing split
            frc_s_i = fit(r_i, frc_i, 'smoothingSpline');
            
            % calculate one/half bit line based on the number of points
            if strcmp(resThreshMethod, 'one-bit')
                thresh_points = (sqrt(npoints(:, i)) + 2 * sqrt(2) + 2) ./ (3 * sqrt(npoints(:, i)) + 2 * sqrt(2));
                thresh_lim = 1 / 3;
            else
                thresh_points = (sqrt(npoints(:, i)) * 0.2071 + 1.9102) ./ (sqrt(npoints(:, i)) * 1.2071 +0.9102);   
                thresh_lim = 0.2071 / 1.2071;
            end
            
            % fit the threshold line to decide the intersection
            thresh_s_i = fit([r_i(2 : end); 2], [thresh_points(2 : end); thresh_lim], 'smoothingSpline');
            
            rs_i = r_i(1) : 0.0001 : 1;
            diff_s_i = frc_s_i(rs_i) - thresh_s_i(rs_i);
            ind = find(diff_s_i <= 0, 1, 'first');
            if isempty(ind)
                rs_i = r_i(1) : 0.0001 : sqrt(2);                
                frc_s_i = fit(r_i, frc_i, 'poly2');
                diff_s_i = frc_s_i(rs_i) - thresh_s_i(rs_i);
                ind = find(diff_s_i <= 0, 1, 'first');
            end
            if isempty(ind)
                points(i) = sqrt(2);
                % points(i) = nan;
            else 
                points(i) = rs_i(ind);
            end
        end  
end

psz = xyPixelSize + (dz - xyPixelSize)  * abs(cos(thetas));

res_mat = 2 * psz' ./ points';
dtheta = diff(thetas);
dtheta = dtheta(1);
res_theta = 0 : dtheta: 2 * pi;

% convert the output such that 0 is for x, and pi/2 for z
cind = (numel(res_mat)) / 2 + 1;
res_out = [res_mat(cind : end), res_mat(1:end), res_mat(1 : cind)]';

res.comp_angles = thetas; 
res.points = points; 
res.thetas = res_theta';
res.res = res_out;


% test cutoff correction (from miplib), we may need our own calibration for
% our images.
if false
    res_1_mat = res_mat;

    params = [0.95988146, 0.97979108, 13.90441896, 0.55146136];

    for i = 1 : numel(res_mat)
        point = points(i);

        cut_off_correction = params(1) * exp(params(3) * (point - params(2))) + params(4);
        res_1_mat(i) = res_mat(i) / cut_off_correction;
    end

    figure, polarplot(res_theta,[res_1_mat(cind : end), res_1_mat(1:end), res_1_mat(1 : cind)])
    res_out_1 = [res_1_mat(cind : end), res_1_mat(1:end), res_1_mat(1 : cind)]';
    
    res.res = res_out_1;
end

end

