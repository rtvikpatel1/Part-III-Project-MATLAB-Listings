clc; clear; close all;

%% === PARAMETERS ===
impingement_angle = 30;             % degrees
Q_LPM = 2;                          % Flow rate in LPM
Q = Q_LPM / 1000 / 60;              % Flow rate in m³/s
nu = 1.0016e-6;                     % Kinematic viscosity (m²/s)
rho = 1000;                         % Density (kg/m³)
mu = nu * rho;                      % Dynamic viscosity (Pa·s)
gamma = 72.8e-3;                    % Surface tension (N/m)
g = 9.81;                           % Gravity (m/s²)
h0 = 0.1e-3;
U_0 = Q/(pi*(1.5e-3)^2);
                    % Initial film height (m)

theta_degrees = 0:1:360;
r_span = [1e-3, 175e-3];            % Radial span [m]
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

% Elliptical impingement radius function (1.5 mm base radius)
r_e = @(theta) 1.5e-3 .* (sind(90 - impingement_angle) ./ ...
    (1 - cosd(theta) .* cosd(90 - impingement_angle)));

% Angular flow rate function q(θ)
Q_theta = @(theta) Q .* (sind(90 - impingement_angle)^3 ./ ...
    (1 - cosd(theta) .* cosd(90 - impingement_angle)).^2);

A = 0.2964;

% Effective Reynolds number and transition radius

%r_b_3 = @(theta) r_e(theta).^4 .* sind(90-impingement_angle).^2 * (U_0 * nu) / (2.12)^2 ;
%B = @(theta) 1.9 * mu /(rho*(U_0*r_e(theta).^2 * sind(90-impingement_angle)).^2);
%r_crit_theta = @(theta) (1/(U_0*B(theta)) - r_b_3(theta)/2)^(1/3);


Re_eff_theta = @(theta) (1/nu) * (U_0) / (2*r_e(theta));
r_crit_theta = @(theta) A .* 3e-3 .* (Re_eff_theta(theta)).^(1/3);

%% === MODEL CONSTANTS ===
% Laminar
C1_1 = 0.6137; C2_1 = 0.4755; f_prime_0_1 = 1.402;

% Turbulent
C1_2 = 0.875; C2_2 = 0.778;

%% === STORAGE ARRAYS ===
radii_laminar = NaN(size(theta_degrees));
radii_turbulent = NaN(size(theta_degrees));
radii_hybrid = NaN(size(theta_degrees));

%% === MAIN LOOP ===
for j = 1:numel(theta_degrees)
    theta = theta_degrees(j);
    Q_star = Q_theta(theta);
    r_crit = r_crit_theta(theta);

    % --- Define ODEs ---
    F1 = @(r, h) ...
        (mu * f_prime_0_1 * Q_star ./ (2 * pi * C1_1 .* r .* h.^2) ...
        - C2_1 * rho * (Q_star / (2 * pi * C1_1))^2 ./ (r.^3 .* h)) ...
        ./ (C2_1 * rho * (Q_star / (2 * pi * C1_1))^2 ./ (r.^2 .* h.^2) ...
        - gamma ./ h - 0.5 * rho * g .* h);

    F2 = @(r, h) ...
        (rho * ((Q_star ./ (2 * pi * C1_2 .* r .* h)).^2) .* ...
        0.225 .* ((mu * 2 * pi * C1_2 .* r) ./ (rho * Q)).^(1/4) ...
        - C2_2 * rho * (Q_star / (2 * pi * C1_2))^2 ./ (r.^3 .* h)) ...
        ./ (C2_2 * rho * (Q_star / (2 * pi * C1_2))^2 ./ (r.^2 .* h.^2) ...
        - gamma ./ h - 0.5 * rho * g .* h);

    %% === Pure Laminar Model ===
    try
        [rL, hL] = ode15s(F1, r_span, h0, options);
        if numel(rL) > 5
            dh = diff(hL) ./ diff(rL);
            idx = find(abs(dh) > 5000, 1);
            if isempty(idx)
                radii_laminar(j) = rL(end);
            else
                radii_laminar(j) = rL(idx);
            end
        end
    catch
        warning("Laminar model failed at θ = %g°", theta);
    end

    %% === Pure Turbulent Model ===
    try
        [rT, hT] = ode15s(F2, r_span, h0, options);
        if numel(rT) > 5
            dh = diff(hT) ./ diff(rT);
            idx = find(abs(dh) > 5000, 1);
            if isempty(idx)
                radii_turbulent(j) = rT(end);
            else
                radii_turbulent(j) = rT(idx);
            end
        end
    catch    
        warning("Turbulent model failed at θ = %g°", theta);
    end

    %% === Hybrid Model ===
    event_transition = @(r, h) deal_transition_rcrit(r, r_crit);
    options_with_event = odeset(options, 'Events', event_transition);

    try
        [r1, h1] = ode15s(F1, r_span, h0, options_with_event);
        if isempty(r1), continue; end
        r_trans = r1(end);
        h_transition = h1(end);

        [r2, h2] = ode15s(F2, [r_trans, r_span(2)], h_transition, options);

        r_full = [r1; r2];
        h_full = [h1; h2];

        if numel(r_full) > 5
            dh = diff(h_full) ./ diff(r_full);
            idx = find(abs(dh) > 5000, 1);
            if isempty(idx)
                radii_hybrid(j) = r_full(end);
            else
                radii_hybrid(j) = r_full(idx);
            end
        end
    catch
        warning("Hybrid model failed at θ = %g°", theta);
    end
end

%% === PLOT RESULTS ===
figure('Color', 'w');
polarplot(deg2rad(theta_degrees), radii_laminar*1000, 'b-', 'LineWidth', 2); hold on;
polarplot(deg2rad(theta_degrees), radii_turbulent*1000, 'r--', 'LineWidth', 2);
polarplot(deg2rad(theta_degrees), radii_hybrid*1000, 'k:', 'LineWidth', 2);

% Formatting
ax = gca;
ax.ThetaTick = 0:30:360;
ax.ThetaLim = [0 360];
ax.RLim = [0, max([radii_laminar, radii_turbulent, radii_hybrid]*1000) * 1.1];
ax.FontSize = 14;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';

legend({'Laminar', 'Turbulent', 'Hybrid'}, 'Location', 'southoutside', 'FontSize', 12);
title('Hydraulic Jump Radius Comparison (0°–360°)', 'FontSize', 16);

%% === TRANSITION FUNCTION ===
function [value, isterminal, direction] = deal_transition_rcrit(r, r_crit)
    value = r - r_crit;     % Event when radius reaches critical value
    isterminal = 1;         % Stop integration
    direction = 1;          % Trigger when crossing 0 from below
end
