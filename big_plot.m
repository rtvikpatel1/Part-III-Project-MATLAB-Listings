clc; clear; close all;

%% === USER CONFIGURATION ===
dfi_main = 'your_main_image.dfi';
dfi_bg   = 'your_background.dfi';
csv_out  = 'jump_radius_data.csv';  % final spreadsheet

scale_mm_per_pixel = 0.1151;
absorption_coeff   = 2.1988;

impingement_angle = 20;  % degrees
Q_LPM = 1.19;
nu = 20.76e-6;
rho = 1160;
gamma = 67e-3;
g = 9.81;

angles_deg = 0:30:180;

%% === LOAD IMAGE AND COMPUTE HEIGHT MAP ===
main_img = double(dfi2mat(dfi_main).image);
bg_img   = double(dfi2mat(dfi_bg).image);
div_img  = main_img ./ bg_img;
div_img(div_img <= 0 | isnan(div_img)) = eps;
height_map = -absorption_coeff * log(div_img);

%% === SELECT CENTER ===
figure;
imshow(height_map, stretchlim(height_map, [0.01 0.99]), 'InitialMagnification', 'fit');
title('Click center of jump'); colormap gray;
[cx, cy] = ginput(1);
close;

%% === COLLECT MEASURED RADII ===
r_max_px = round(100 / scale_mm_per_pixel);
jump_radii_px = nan(size(angles_deg));

figure;
imshow(height_map, stretchlim(height_map, [0.01 0.99]), 'InitialMagnification', 'fit');
colormap gray;
title('Click jump radius on each blue line'); hold on;
plot(cx, cy, 'ro');

theta_rad = deg2rad(angles_deg);
for i = 1:length(theta_rad)
    x_end = cx + r_max_px * cos(theta_rad(i));
    y_end = cy - r_max_px * sin(theta_rad(i));
    line([cx x_end], [cy y_end], 'Color', 'b');
end

for i = 1:length(theta_rad)
    [x_jump, y_jump] = ginput(1);
    dx = x_jump - cx;
    dy = y_jump - cy;
    jump_radii_px(i) = sqrt(dx^2 + dy^2);
end

measured_mm = jump_radii_px * scale_mm_per_pixel;

%% === THEORETICAL RADIUS: BHAGAT 2020 ===
Q = Q_LPM / 1000 / 60;
c = cosd(90 - impingement_angle);
s = sind(90 - impingement_angle);

q = @(theta) (Q / ((12e4) * pi)) * (s^3 ./ max((1 + c * cosd(180 - theta)).^2, 1e-6));
Q_star = @(theta) gamma^2 ./ (2*pi*q(theta) .* rho^2 .* nu .* g);
r_mixed_model = @(theta) 0.2705 .* (2*pi)^(3/4) .* q(theta).^(3/4) .* ...
    nu^(-1/4) .* gamma^(-1/4) .* rho^(1/4) .* ...
    (sqrt(Q_star(theta).^2 + 2.*Q_star(theta)) - Q_star(theta)).^0.25;

theory_mm = r_mixed_model(angles_deg) * 1000;

%% === APPEND TO CSV ===
n = length(angles_deg);
filename_tag = repmat({dfi_main}, n, 1);

T = table(filename_tag, angles_deg(:), measured_mm(:), theory_mm(:), ...
    'VariableNames', {'Filename', 'Angle_deg', 'Measured_mm', 'Theory_mm'});

if isfile(csv_out)
    T_existing = readtable(csv_out);
    T_combined = [T_existing; T];
    writetable(T_combined, csv_out);
else
    writetable(T, csv_out);
end

disp('✅ Data saved to spreadsheet:');
disp(csv_out);
