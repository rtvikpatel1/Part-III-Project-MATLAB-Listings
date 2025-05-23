clc; clear; close all;

%% === USER SETTINGS ===
dfi_main = 'rota_7_avg_15deg_3mm_nozzle.dfi';
dfi_bg = 'background_avg_15deg_3mm_nozzle.dfi';
scale_mm_per_pixel = 0.11828;

rota_value = 7;
impingement_angle = 15;
Q = 0.5 + (rota_value - 2)*0.3;    % in LPM
nozzle_radius = 1.5e-3;
nu = 1.0016e-6;
rho = 1e3;
gamma = 72.8e-3;
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






% --- Physical Parameters ---
nu = 1.0016e-6;
rho = 1000;
mu = nu * rho;
gamma = 72.8e-3;
g = 9.81;
C1 = 0.6137;
C2 = 0.4755;
f_prime_0 = 1.402;
Q_LPM = 2;
Q = Q_LPM / 1000 / 60;

% --- Setup ---
theta_degrees = 0:1:360;                     % Angular steps (adjustable)
r_span = [15e-3, 150e-3];                    % Radius span [m]
r_vals = linspace(r_span(1), r_span(2), 300);
h0 = 0.1e-3;                                 % Initial height
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

% Angular flow rate function
Q_theta = @(theta) Q * ((sind(90 - impingement_angle)^3) ./ ((1 - cosd(theta) * cosd(90 - impingement_angle)).^2));

% Preallocate
num_angles = numel(theta_degrees);
h_grid = NaN(num_angles, numel(r_vals));
singular_radii = NaN(1, num_angles);

% --- Solve ODE for each angle ---
for j = 1:num_angles
    theta = theta_degrees(j);
    Q_star = Q_theta(theta);

    F = @(r, h) ...
        (mu * f_prime_0 * Q_star ./ (2 * pi * C1 .* r .* h.^2) ...
        - C2 * rho * (Q_star / (2 * pi * C1))^2 ./ (r.^3 .* h)) ...
        ./ (C2 * rho * (Q_star / (2 * pi * C1))^2 ./ (r.^2 .* h.^2) ...
        - gamma ./ h - 0.5 * rho * g .* h);

    try
        [r_sol, h_sol] = ode15s(F, r_span, h0, options);

        if numel(r_sol) > 5
            dh = diff(h_sol) ./ diff(r_sol);
            slope_threshold = 5000;
            unstable_idx = find(abs(dh) > slope_threshold, 1);

            if ~isempty(unstable_idx)
                r_cut = r_sol(unstable_idx);
            else
                r_cut = r_sol(end) - 2e-3;
            end

            h_interp = interp1(r_sol, h_sol, r_vals, 'linear', NaN);
            h_interp(r_vals > r_cut) = NaN;
            h_grid(j, :) = h_interp;
            singular_radii(j) = r_cut;
        end
    catch
        warning("ODE failed at θ = %g°", theta);
    end
end





%% === POLAR PLOT 0°–180° WITH STYLING ===

%% === POLAR PLOT 0°–180° WITH FULL FORMATTING ===
figure('Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.6 0.8]);  % large figure with white background

% Plot curves
p1 = polarplot(deg2rad(theta_deg), jump_radii_mm, 'go', 'LineWidth', 2.5, 'MarkerSize', 12); hold on;
p2 = polarplot(deg2rad(theta_degrees), singular_radii * 1000, 'LineWidth', 2);

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
legend({'Experiment', 'Bhagat et al. (2022)'}, ...
    'Interpreter', 'latex', 'FontSize', 25, 'Location', 'southoutside');


% Add border to figure
annotation('rectangle', [0 0 1 1], 'Color', 'k', 'LineWidth', 2);  % border frame



%% === Auto-generate filename and save ===
% Create a filename string using parameters
filename_base = sprintf('jump_location_integrated_theoryrota_%d_%d_deg_3mm_nozzle', ...
    rota_value, impingement_angle);

% Save as PNG
saveas(gcf, [filename_base '.png']);



