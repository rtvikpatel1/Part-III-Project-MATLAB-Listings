clc; clear; close all;

%% === USER SETTINGS ===
dfi_main = 'avg_rota_highest_20_deg_glycerol.dfi';
dfi_bg = 'avg_background_20_deg_glycerol.dfi';
scale_mm_per_pixel = 0.1151;

rota_value = 'highest';
impingement_angle = 20;
Q = 1.19;    % in LPM
nozzle_radius = 1.5e-3;
nu = 20.76e-6;
rho = 1160;
gamma = 67e-3;
g = 9.81;



angle_step = 15;                     % Angular resolution
r_max_mm = 100;                     % Max radius for guide lines
r_max_px = round(r_max_mm / scale_mm_per_pixel);


%% === LOAD & COMPUTE HEIGHT MAP ===
main_img = double(dfi2mat(dfi_main).image);
bg_img = double(dfi2mat(dfi_bg).image);
div_img = main_img ./ bg_img;
div_img(div_img <= 0) = NaN;
height_map = -2.1988 * log(div_img);
[rows, cols] = size(height_map);

%% === SELECT CENTER ===
figure;
contrast_limits = stretchlim(height_map, [0.01 0.99]);  % Adjust contrast using 1st–99th percentile
imshow(height_map, contrast_limits, 'InitialMagnification', 'fit');
colormap gray;
title('Click on the center of the jump');
[cx, cy] = ginput(1);
center = [cx, cy];

% Re-plot and zoom in for jump radius selection
figure;
imshow(height_map, contrast_limits, 'InitialMagnification', 'fit');
colormap gray;
title('Click the jump radius on each blue line');
hold on;
plot(cx, cy, 'ro');

% Zoomed-in view around center (adjust this size as needed)
zoom_half_width = 150;  % pixels
xlim([cx - zoom_half_width, cx + zoom_half_width]);
ylim([cy - zoom_half_width, cy + zoom_half_width]);



%% === SAMPLE RADIUS FROM 0° TO 180° ===
theta_deg = 0:angle_step:180;           % For labeling and plots
num_angles = length(theta_deg);
jump_radii_px = NaN(size(theta_deg));

figure;
imshow(height_map, contrast_limits, 'InitialMagnification', 'fit');
colormap gray;
title('Click the jump radius on each blue line');
hold on;
plot(cx, cy, 'ro');

for i = 1:num_angles
    theta_rad = deg2rad(theta_deg(i));
    x_end = cx + r_max_px * cos(theta_rad);
    y_end = cy - r_max_px * sin(theta_rad);

    line([cx x_end], [cy y_end], 'Color', 'b');
    drawnow;

    [x_jump, y_jump] = ginput(1);
    dx = x_jump - cx;
    dy = y_jump - cy;
    r_px = sqrt(dx^2 + dy^2);
    jump_radii_px(i) = r_px;
end



%% === CONVERT TO MM AND EVALUATE THEORY ===
jump_radii_mm = jump_radii_px * scale_mm_per_pixel;

% Define your theoretical function here
c = cosd(90 - impingement_angle);  % cos(angle from jet to plate)
s = sind(90 - impingement_angle);

q = @(theta) (Q/((12e4)*pi)) * (s^3 ./((1 + c*cosd(180-theta)).^2));

theta_smooth = 0:1:180;  % finer resolution



%% === GRAVITY-DOMINATED THEORETICAL MODEL ===
% Based on Bohr et al. (1993) style model
r_gravity_model = @(theta) 0.73 .* q(theta).^(5/8) .* nu^(-3/8) .* g^-(1/8);

r_gravity_mm = r_gravity_model(theta_smooth) .* 1000;



%% === SURFACE TENSION-DOMINATED THEORETICAL MODEL ===
r_surface_model = @(theta) 0.2705 .* (2*pi)^(3/4) .* q(theta).^(3/4) .* nu^(-1/4) .* gamma^(-1/4) .* rho^(1/4);

r_surface_mm = r_surface_model(theta_smooth) .* 1000;


%% === GRAVITY+ ST DOMINATED THEORETICAL MODEL ===
Q_star = @(theta) gamma^2 ./ (2*pi*q(theta) .* rho^2 .* nu .* g);
r_mixed_model = @(theta) 0.2705 .* (2*pi)^(3/4) .* q(theta).^(3/4) .* nu^(-1/4) .* gamma^(-1/4) .* rho^(1/4) .* (sqrt(Q_star(theta).^2 + 2.*Q_star(theta))-Q_star(theta)).^0.25;

r_mixed_mm = r_mixed_model(theta_smooth) .* 1000;





%% === POLAR PLOT 0°–180° WITH STYLING ===

%% === POLAR PLOT 0°–180° WITH FULL FORMATTING ===
figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.6 0.8]);  % large figure with white background

% Plot curves
p1 = polarplot(deg2rad(theta_deg), jump_radii_mm, 'go', 'LineWidth', 2.5, 'MarkerSize', 12); hold on;
p2 = polarplot(deg2rad(theta_smooth), r_gravity_mm, 'k-', 'LineWidth', 2);
p3 = polarplot(deg2rad(theta_smooth), r_surface_mm, 'b-', 'LineWidth', 2);
p4 = polarplot(deg2rad(theta_smooth), r_mixed_mm, 'r-', 'LineWidth', 2);

% Axis styling
ax = gca;
ax.ThetaTick = 0:30:180;
ax.ThetaTickLabel = {'$0^\circ$','$30^\circ$','$60^\circ$','$90^\circ$', ...
                     '$120^\circ$','$150^\circ$','$180^\circ$'};
ax.ThetaColor = 'k';
ax.RColor = 'k';
ax.LineWidth = 1.5;
ax.FontSize = 30;
ax.TickLabelInterpreter = 'latex';
thetalim([0 180]);

% Force axes to update
drawnow;

% === Adjust polar axes position to give legend room ===
ax = gca;
ax.Position = [0.15 0.25 0.65 0.65];  % [left bottom width height]

% === Draw border around the polar plot (larger than ax) ===
borderPadding = 0.08;  % how much larger than ax.Position
borderPos = [
    ax.Position(1) - borderPadding, ...
    ax.Position(2) - borderPadding, ...
    ax.Position(3) + 2*borderPadding, ...
    ax.Position(4) + 2*borderPadding];
annotation('rectangle', borderPos, 'Color', 'k', 'LineWidth', 2);

% === Add Radius (mm) label just inside bottom border ===
annotation('textbox', ...
    [0.63, 0.35, 0.20, 0.05], ...
    'String', '$\mathrm{Radius\ (mm)}$', ...
    'Interpreter', 'latex', ...
    'FontSize', 30, ...
    'HorizontalAlignment', 'center', ...
    'EdgeColor', 'none');


% Radius limit
rlim([0 max([jump_radii_mm, r_gravity_mm, r_surface_mm, r_mixed_mm]) * 1.1]);

% Legend
legend({'Experiment', 'Bohr et al. (1993)', 'Bhagat et al. (2018)', 'Bhagat et al. (2020)'}, ...
    'Interpreter', 'latex', 'FontSize', 25, 'Location', 'southoutside');


% Add border to figure
annotation('rectangle', [0 0 1 1], 'Color', 'k', 'LineWidth', 2);  % border frame



%% === Auto-generate filename and save ===
% Create a filename string using parameters
filename_base = sprintf('jump_location_rota_%d_%d_deg_3mm_nozzle', ...
    rota_value, impingement_angle);

% Save as PNG
saveas(gcf, [filename_base '.png']);



