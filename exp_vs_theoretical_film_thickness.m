clc;

%% Parameters
rota_value = 7;
impingement_angle = 20;

scale_mm_per_pixel = 0.11560;
absorption_coeff = 2.5144;
line_length_pixels = 3000;
angles_deg = [0, 30, 60, 90, 120];
colors = {'k', 'b', 'r', 'g', 'm'};

%% Load images
background = double(dfi2mat('background_avg_20deg_3mm_nozzle.dfi').image);
main = double(dfi2mat('rota_7_avg_20deg_3mm_nozzle.dfi').image);

divided = main ./ background;
divided(divided <= 0) = NaN;
height_map = -absorption_coeff * log(divided);  % mm

%% Ask user for center point
figure; imshow(height_map, []);
title('Click to set the center point');
[center_x, center_y] = ginput(1);
hold on; plot(center_x, center_y, 'ro');

%% Initialize figure
figure('Units', 'centimeters', 'Position', [2, 2, 28, 19]); hold on;

legend_entries = {};
legend_handles = [];

%% --- Plot Experimental Profiles ---
for i = 1:length(angles_deg)
    angle_rad = deg2rad(angles_deg(i));
    x1 = center_x + line_length_pixels * cos(angle_rad);
    y1 = center_y - line_length_pixels * sin(angle_rad);

    num_points = round(line_length_pixels);
    x_line = linspace(center_x, x1, num_points);
    y_line = linspace(center_y, y1, num_points);

    profile = zeros(num_points, 1);
    for j = 1:num_points
        x_pix = round(x_line(j));
        y_pix = round(y_line(j));
        if x_pix >= 1 && x_pix <= size(height_map, 2) && y_pix >= 1 && y_pix <= size(height_map, 1)
            profile(j) = height_map(y_pix, x_pix);
        else
            profile(j) = NaN;
        end
    end

    radius = (0:num_points-1) * scale_mm_per_pixel;

    % Mask out certain radius intervals
    if angles_deg(i) == 120
        profile(radius >= 9.5 & radius <= 13 ) = NaN;
        profile(radius >= 22 & radius <= 25 ) = NaN;
        profile(radius >= 65 ) = NaN;
    elseif angles_deg(i) == 30
        profile(radius >= 20 & radius <= 26 ) = NaN;
    end

    h_exp = plot(radius, profile, ':', 'LineWidth', 2, 'Color', colors{i}); % Dotted line

    legend_handles(end+1) = h_exp;
    legend_entries{end+1} = sprintf('$%d^\\circ$ Exp', angles_deg(i));
end

%% --- Plot Theoretical Profiles ---
% Physical parameters
nu = 1.0016e-6;
rho = 1000;
mu = nu * rho;
gamma = 72.8e-3;
g = 9.81;
C1 = 0.6137;
C2 = 0.4755;
f_prime_0 = 1.402;
Q_LPM = 2.0;
imp_ang = 20;  % Use consistent impingement angle for theory
Q = Q_LPM / 1000 / 60;
h0 = 0.05e-3;
r_span = [10e-3, 170e-3];

Q_theta = @(theta) Q * ((sind(90-imp_ang)^3) / ((1 - cosd(theta) * cosd(90-imp_ang))^2));
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

for j = 1:length(angles_deg)
    theta = angles_deg(j);
    Q_star = Q_theta(theta);

    F = @(r, h) ...
        (mu * f_prime_0 * Q_star ./ (2*pi*C1 .* r .* h.^2) ...
        - C2 * rho * (Q_star / (2*pi*C1))^2 ./ (r.^3 .* h)) ...
        ./ (C2 * rho * (Q_star / (2*pi*C1))^2 ./ (r.^2 .* h.^2) ...
        - gamma ./ h - 0.5 * rho * g .* h);

    try
        [r_sol, h_sol] = ode15s(F, r_span, h0, options);
        if length(r_sol) > 1
            h_th = plot(r_sol*1000, h_sol*1000, '-', 'Color', colors{j}, 'LineWidth', 2); % Solid line
            legend_handles(end+1) = h_th;
            legend_entries{end+1} = sprintf('$%d^\\circ$ Theory', angles_deg(j));
        end
    catch
        warning('Failed for angle %d°', theta);
    end
end

%% Final formatting
xlabel('Radius (mm)', 'Interpreter', 'latex', 'FontSize', 30);
ylabel('Film Thickness (mm)', 'Interpreter', 'latex', 'FontSize', 30);


grid on; box on;
set(gca, 'FontSize', 25, 'TickLabelInterpreter', 'latex');
xlim([0 89]);
ylim([0 4]);

%% Save
filename_base = sprintf('exp_vs_theory_profiles_rota_%d_%d_deg_3mm_nozzle', rota_value, impingement_angle);
saveas(gcf, [filename_base '.png']);
