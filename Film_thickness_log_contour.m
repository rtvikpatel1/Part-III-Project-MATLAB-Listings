clc
% clear all
% close all

% Define a list of colors to cycle through
colors = {'k', 'b', 'r', 'g'}; % k: black, b: blue, r: red, g: green

% Initialize color index if it does not exist
if ~exist('color_index', 'var')
    color_index = 1;
end
clc;
% clear all;
% close all;

%% Settings
scale_mm_per_pixel = 0.1151; % Modify as per your image scale
absorption_coeff = 2.6645;
rota_value = 'highest';
impingement_angle = 0;


%% Step 1: Load images
background_image_dfi = dfi2mat('avg_background_0_deg_glycerol.dfi');
main_image_dfi = dfi2mat('avg_rota_highest_0_deg_glycerol.dfi');
background_image = double(background_image_dfi.image);
main_image = double(main_image_dfi.image);

%% Step 2: Divide images and compute film height
divided_image = main_image ./ background_image;

% Avoid log of 0 or negative values
divided_image(divided_image <= 0) = NaN;

% Film thickness (height) computation
height_map = -absorption_coeff * log(divided_image); % In mm

%% Step 3: Rescale axes for real-world dimensions
[img_height, img_width] = size(height_map);
x_mm = (0:img_width-1) * scale_mm_per_pixel;
y_mm = (0:img_height-1) * scale_mm_per_pixel;

%% Step 4: Apply log10 and create contour with log-scaled appearance
height_map(height_map <= 0.01) = 0.01; % Avoid log(0) and negatives
height_map(height_map > 5) = 5;        % Clamp above 5 mm

log_height_map = log10(height_map);    % Take log10 for contour

% Define log-spaced contour levels (log10 values)
log_levels = linspace(log10(0.01), log10(5), 100);

figure('Units', 'centimeters', 'Position', [2, 2, 28, 19]);
[~, h] = contourf(x_mm, y_mm, log_height_map, log_levels);
set(h, 'LineColor', 'none');
colormap jet;

% Create colorbar and manually set ticks for actual values
c = colorbar;
caxis([log10(0.1), log10(5)]);

% Set ticks at meaningful powers of 10
tick_values = log10([0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5]);
% Generate LaTeX-formatted tick labels
tick_labels = arrayfun(@(x) sprintf('$%.2g$', 10^x), tick_values, 'UniformOutput', false);

% Apply to colorbar
set(c, 'Ticks', tick_values, 'TickLabels', tick_labels, 'TickLabelInterpreter', 'latex');
c.Label.String = 'Film Thickness (mm)';
c.Label.Interpreter = 'latex';
c.Label.FontSize = 25;

xlabel('X (mm)', 'Interpreter', 'latex', 'FontSize', 30);
ylabel('Y (mm)', 'Interpreter', 'latex', 'FontSize', 30);
set(gca, 'FontSize', 25, 'TickLabelInterpreter', 'latex');
axis equal tight;
box on;

%% === Auto-generate filename and save ===
% Create a filename string using parameters
filename_base = sprintf('contour_rota_%d_%d_deg_glycerol', ...
    rota_value, impingement_angle);

% Save as PNG
saveas(gcf, [filename_base '.png']);

