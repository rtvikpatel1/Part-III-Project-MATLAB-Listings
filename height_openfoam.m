clc; clear; close all;

%% === USER PARAMETERS ===
dfi_main = 'avg_rota_10_20_deg_glycerol.dfi';
dfi_bg   = 'avg_background_20_deg_glycerol.dfi';
foam_csv = ['20_deg_before.csv'];

absorption_coeff = 2.6645;
scale_mm_per_pixel = 0.11610;
angles_deg = [0, 30, 60, 90, 120];
colors = {'k', 'b', 'r', 'g', 'm'};

r_max_mm = 100;
r_ignore_mm = 2;
height_threshold = 1.1;

%% === LOAD EXPERIMENTAL HEIGHT MAP ===
main = double(dfi2mat(dfi_main).image);
bg = double(dfi2mat(dfi_bg).image);
divided = main ./ bg;
divided(divided <= 0 | isnan(divided)) = eps;
height_map = -absorption_coeff * log(divided); % mm

%% === LOAD OPENFOAM CONTOUR DATA ===
foam = readtable(foam_csv);
Xf = foam.Points_0 * 1000;
Yf = foam.Points_1 * 1000;
Zf = foam.Points_2 * 1000;
F = scatteredInterpolant(Xf, Yf, Zf, 'natural', 'none');

% Meshgrid for visualization & interpolation
xq = linspace(min(Xf), max(Xf), 300);
yq = linspace(min(Yf), max(Yf), 300);
[Xq, Yq] = meshgrid(xq, yq);
Zq = F(Xq, Yq);

%% === APPROVAL-BASED PROFILE COMPARISON ===
accepted_r_exp = {};
accepted_z_exp = {};
accepted_r_foam = {};
accepted_z_foam = {};
accepted_colors = {};

for i = 1:length(angles_deg)
    angle = angles_deg(i);
    angle_rad = deg2rad(angle);
    approved = false;

    while ~approved
        %% Select center on experimental image
        figure(1); clf;
        imshow(height_map, []); colormap gray;
        title(sprintf('Click EXPERIMENTAL center for %d°', angle));
        [cx_exp, cy_exp] = ginput(1);
        hold on; plot(cx_exp, cy_exp, 'ro');

        %% Select center on OpenFOAM surface
        figure(2); clf;
        surf(Xq, Yq, Zq, 'EdgeColor', 'none'); view(2); axis equal;
        colormap turbo; colorbar;
        title(sprintf('Click FOAM center for %d°', angle));
        [cx_foam, cy_foam] = ginput(1);
        close(2);

        %% === Extract experimental profile ===
        r_pixels = round(r_max_mm / scale_mm_per_pixel);
        x_end = cx_exp + r_pixels * cos(angle_rad);
        y_end = cy_exp - r_pixels * sin(angle_rad);  % y flipped for images

        x_line = linspace(cx_exp, x_end, r_pixels);
        y_line = linspace(cy_exp, y_end, r_pixels);
        z_exp = nan(1, r_pixels);

        for j = 1:r_pixels
            x = round(x_line(j));
            y = round(y_line(j));
            if x >= 1 && x <= size(height_map,2) && y >= 1 && y <= size(height_map,1)
                z_exp(j) = height_map(y,x);
            end
        end

        r_exp = (0:r_pixels-1) * scale_mm_per_pixel;

        %% === Extract FOAM profile ===
        r_foam_vec = linspace(0, r_max_mm, 1000);
        x_foam_line = cx_foam + r_foam_vec * cos(angle_rad);
        y_foam_line = cy_foam + r_foam_vec * sin(angle_rad);
        z_foam_line = F(x_foam_line, y_foam_line);

        % Smooth and threshold
        r_valid = r_foam_vec >= r_ignore_mm;
        r_valid_steps = r_foam_vec(r_valid);
        z_valid = z_foam_line(r_valid);

        %idx_jump = find(z_valid > height_threshold, 1, 'first');
        %if ~isempty(idx_jump)
            %r_cutoff = r_valid_steps(idx_jump);
        %else
            %r_cutoff = max(r_foam_vec);
        %end

        %mask = r_foam_vec <= r_cutoff;
        r_foam = r_foam_vec;
        z_foam = smoothdata(z_foam_line, 'movmean', 15);


        %% === Preview plot ===
        figure(3); clf; hold on;
        plot(r_exp, z_exp, '-', 'Color', colors{i}, 'LineWidth', 2);
        plot(r_foam, z_foam, '--', 'Color', colors{i}, 'LineWidth', 2);
        xlabel('Radius (mm)', 'Interpreter', 'latex');
        ylabel('Film Thickness (mm)', 'Interpreter', 'latex');
        title(sprintf('%d° Slice Comparison', angle), 'Interpreter', 'latex');
        grid on; xlim([0 35]); ylim([0 5]);

        %% === Approve? ===
        resp = input('Approve this slice? [y/n]: ', 's');
        if strcmpi(resp, 'y')
            accepted_r_exp{end+1} = r_exp;
            accepted_z_exp{end+1} = z_exp;
            accepted_r_foam{end+1} = r_foam;
            accepted_z_foam{end+1} = z_foam;
            accepted_colors{end+1} = colors{i};
            approved = true;
        end
    end
end

%% === 3D Plot: All Approved Slices ===
figure('Name', '3D Profile Comparison', 'Color', 'w');
hold on;

for i = 1:length(accepted_r_exp)
    % Filter experimental to r <= 75
    mask_exp = accepted_r_exp{i} <= 75;
    r_exp = accepted_r_exp{i}(mask_exp);
    z_exp = accepted_z_exp{i}(mask_exp);
    y_exp = angles_deg(i) * ones(size(r_exp));

    % Filter OpenFOAM to r <= 75
    r_foam = accepted_r_foam{i};
    z_foam = accepted_z_foam{i};
    mask_foam = r_foam <= 75;
    r_foam = r_foam(mask_foam);
    z_foam = z_foam(mask_foam);
    y_foam = angles_deg(i) * ones(size(r_foam));

    % Plot both
    plot3(r_exp, y_exp, z_exp, '-', 'Color', colors{i}, 'LineWidth', 2);
    plot3(r_foam, y_foam, z_foam, '--', 'Color', colors{i}, 'LineWidth', 2);
end

% Axis labels
xlabel('Radius (mm)', 'Interpreter', 'latex');
zlabel('Film Thickness (mm)', 'Interpreter', 'latex');
% No y-label
ylabel('');

% Grid and axes
grid on; box on;
view(3);
set(gca, 'FontSize', 12, 'TickLabelInterpreter', 'latex');

% Y-axis ticks and labels for angles
yticks(angles_deg);
yticklabels({'$0^\circ$', '$30^\circ$', '$60^\circ$', '$90^\circ$', '$120^\circ$'});
set(gca, 'TickLabelInterpreter', 'latex');

% Axis limits
xlim([0 75]);

% Remove title
title('');
saveas(gcf, 'openfoam_3D_height_final.png');

