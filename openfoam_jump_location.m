clc; clear; close all;

%% === USER SETTINGS ===
dfi_main = 'avg_rota_10_0_deg_glycerol.dfi';
dfi_bg   = 'avg_background_0_deg_glycerol.dfi';
scale_mm_per_pixel = 0.1151;
absorption_coeff   = 2.1988;

foam_csv = '0_deg_after.csv';

rota_value         = '10';
impingement_angle  = 0;
Q_LPM              = 1;
nu                 = 20.76e-6;
rho                = 1160;
gamma              = 67e-3;
g                  = 9.81;

angle_step_deg     = 15;
r_max_mm           = 100;
r_ignore_mm        = 5;
height_threshold   = 0.79;

%% === LOAD OPENFOAM SURFACE DATA ===
foam_data = readtable(foam_csv);
X = foam_data.Points_0 * 1000;
Y = foam_data.Points_1 * 1000;
Z = foam_data.Points_2 * 1000;

F = scatteredInterpolant(X, Y, Z, 'natural', 'none');
xq = linspace(min(X), max(X), 300);
yq = linspace(min(Y), max(Y), 300);
[Xq, Yq] = meshgrid(xq, yq);
Zq = F(Xq, Yq);

figure;
surf(Xq, Yq, Zq, 'EdgeColor','none'); view(2); axis equal;
xlabel('$x$ (mm)', 'Interpreter','latex');
ylabel('$y$ (mm)', 'Interpreter','latex');
title('Select Jump Center', 'Interpreter','latex');
colormap turbo; colorbar;
[xc_mm, yc_mm] = ginput(1);
close;

%% === EXTRACT FOAM RADIUS PROFILE ===
angles_deg_foam = 0:5:180;
theta_rad_foam = deg2rad(angles_deg_foam);
r_steps = linspace(0, r_max_mm, 1000);
jump_radii_foam = nan(size(theta_rad_foam));

for i = 1:length(theta_rad_foam)
    theta = theta_rad_foam(i);
    x_line = xc_mm + r_steps .* cos(theta);
    y_line = yc_mm + r_steps .* sin(theta);
    z_line = F(x_line, y_line);

    r_valid = r_steps >= r_ignore_mm;
    r_valid_steps = r_steps(r_valid);
    z_valid = z_line(r_valid);

    idx = find(z_valid > height_threshold, 1, 'first');
    if ~isempty(idx)
        jump_radii_foam(i) = r_valid_steps(idx);
    end
end

valid = ~isnan(jump_radii_foam);
theta_foam_smooth = linspace(min(theta_rad_foam(valid)), max(theta_rad_foam(valid)), 300);
r_foam_smooth = interp1(theta_rad_foam(valid), jump_radii_foam(valid), theta_foam_smooth, 'spline');

%% === LOAD EXPERIMENTAL HEIGHT MAP ===
main_img = double(dfi2mat(dfi_main).image);
bg_img   = double(dfi2mat(dfi_bg).image);
div_img  = main_img ./ bg_img;
div_img(div_img <= 0 | isnan(div_img)) = eps;
height_map = -absorption_coeff * log(div_img);

figure;
contrast_limits = stretchlim(height_map, [0.01 0.99]);
imshow(height_map, contrast_limits, 'InitialMagnification', 'fit');
colormap gray;
title('Click center of jump');
[cx, cy] = ginput(1);
close;

r_max_px = round(r_max_mm / scale_mm_per_pixel);
theta_deg = 0:angle_step_deg:180;
theta_rad = deg2rad(theta_deg);
jump_radii_px = nan(size(theta_deg));

figure;
imshow(height_map, contrast_limits, 'InitialMagnification', 'fit');
colormap gray;
title('Click jump radius on each blue line');
hold on; plot(cx, cy, 'ro');

for i = 1:length(theta_deg)
    x_end = cx + r_max_px * cos(theta_rad(i));
    y_end = cy - r_max_px * sin(theta_rad(i));
    line([cx x_end], [cy y_end], 'Color', 'b');
end

for i = 1:length(theta_deg)
    [x_jump, y_jump] = ginput(1);
    dx = x_jump - cx;
    dy = y_jump - cy;
    jump_radii_px(i) = sqrt(dx^2 + dy^2);
end

jump_radii_mm = jump_radii_px * scale_mm_per_pixel;

%% === THEORY: BHAGAT 2020 MIXED MODEL ===
Q = Q_LPM / 1000 / 60;
c = cosd(90 - impingement_angle);
s = sind(90 - impingement_angle);
q = @(theta) (Q / ((2) * pi)) * (s^3 ./ ((1 + c * cosd(180 - theta)).^2));
theta_smooth_deg = linspace(0, 180, 300);
Q_star = @(theta) gamma^2 ./ (2*pi*q(theta) .* rho^2 .* nu .* g);
r_mixed_model = @(theta) 0.2705 .* (2*pi)^(3/4) .* q(theta).^(3/4) .* ...
    nu^(-1/4) .* gamma^(-1/4) .* rho^(1/4) .* ...
    (sqrt(Q_star(theta).^2 + 2.*Q_star(theta)) - Q_star(theta)).^0.25;
r_mixed_mm = r_mixed_model(theta_smooth_deg) * 1000;

%% === POLAR PLOT ===
figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.6 0.8]);

polarplot(theta_rad, jump_radii_mm, 'go', 'LineWidth', 2.5, 'MarkerSize', 10); hold on;
polarplot(theta_foam_smooth, r_foam_smooth, 'b-', 'LineWidth', 3);
polarplot(deg2rad(theta_smooth_deg), r_mixed_mm, 'r-', 'LineWidth', 2.5);

ax = gca;
ax.ThetaTick = 0:30:180;
ax.ThetaTickLabel = {'$0^\circ$','$30^\circ$','$60^\circ$','$90^\circ$', ...
                     '$120^\circ$','$150^\circ$','$180^\circ$'};
ax.ThetaColor = 'k'; ax.RColor = 'k';
ax.LineWidth = 1.5; ax.FontSize = 30;
ax.TickLabelInterpreter = 'latex';
thetalim([0 180]);
rlim([0 max([jump_radii_mm, r_foam_smooth, r_mixed_mm]) * 1.1]);
ax.Position = [0.15 0.25 0.65 0.65];

annotation('textbox', [0.63, 0.35, 0.20, 0.05], ...
    'String', '$\mathrm{Radius\ (mm)}$', ...
    'Interpreter', 'latex', 'FontSize', 30, ...
    'HorizontalAlignment', 'center', 'EdgeColor', 'none');
annotation('rectangle', [0.07 0.08 0.86 0.86], 'Color', 'k', 'LineWidth', 2);

legend({'Experiment', 'Simulation', 'Bhagat et al. (2020)'}, ...
    'Interpreter', 'latex', 'FontSize', 25, 'Location', 'southoutside');

saveas(gcf, 'jump_radius_combined.png');
