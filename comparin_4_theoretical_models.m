% --- Physical Parameters ---
nu = 1.0016e-6;
rho = 1000;
mu = nu * rho;
gamma = 72.8e-3;
g = 9.81;
C1 = 0.6137;
C2 = 0.4755;
f_prime_0 = 1.402;
rota_value = 7;
Q_LPM = 2;
imp_ang = 45;
Q = Q_LPM / 1000 / 60;

% --- Setup ---
theta_degrees = 0:1:360;                     % Angular steps for ODE
theta_smooth = 0:1:360;                      % Fine angular range for theory curves
r_span = [1e-3, 250e-3];                     % Radius span [m]
r_vals = linspace(r_span(1), r_span(2), 300);
h0 = 0.1e-3;                                 % Initial height
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

% --- Geometry-based Angular Flow Rate Function q(θ) ---
c = cosd(90 - imp_ang);  % Horizontal component
s = sind(90 - imp_ang);  % Vertical component
q = @(theta) (Q .* (s^3 ./ ((1 + c * cosd(180 - theta)).^2)));

% --- Solve ODE for each angle ---
num_angles = numel(theta_degrees);
h_grid = NaN(num_angles, numel(r_vals));
singular_radii = NaN(1, num_angles);

for j = 1:num_angles
    theta = theta_degrees(j);
    Q_star = q(theta);

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

% Convert to mm
singular_radii_mm = singular_radii .* 1000;

%% === THEORETICAL MODELS (GRAVITY, SURFACE, MIXED) ===
q_s = @(theta) (Q/(2*pi) * (s^3 ./((1 + c*cosd(180-theta)).^2)));

% Gravity-dominated
r_gravity_model = @(theta) 0.73 .* q_s(theta).^(5/8) .* nu^(-3/8) .* g^(-1/8);
r_gravity_mm = r_gravity_model(theta_smooth) * 1000;

% Surface tension-dominated
r_surface_model = @(theta) 0.2705 .* (2*pi)^(3/4) .* q_s(theta).^(3/4) .* nu^(-1/4) .* gamma^(-1/4) .* rho^(1/4);
r_surface_mm = r_surface_model(theta_smooth) * 1000;

% Mixed gravity + surface tension model
Q_star_fun = @(theta) gamma^2 ./ (2*pi*q_s(theta) .* rho^2 .* nu .* g);
r_mixed_model = @(theta) 0.2705 .* (2*pi)^(3/4) .* q_s(theta).^(3/4) .* nu^(-1/4) .* gamma^(-1/4) .* rho^(1/4) ...
                     .* (sqrt(Q_star_fun(theta).^2 + 2.*Q_star_fun(theta)) - Q_star_fun(theta)).^0.25;
r_mixed_mm = r_mixed_model(theta_smooth) * 1000;

%% --- 2D Polar Plot of Singularity Radius and Models ---
figure;
polarplot(deg2rad(theta_degrees), singular_radii_mm, 'g-',  'LineWidth', 1); hold on;
polarplot(deg2rad(theta_smooth), r_gravity_mm, 'k-', 'LineWidth', 1);
polarplot(deg2rad(theta_smooth), r_surface_mm, 'b-', 'LineWidth', 1);
polarplot(deg2rad(theta_smooth), r_mixed_mm, 'r-', 'LineWidth', 1);

title('Singularity Radius vs. Angle');
legend('ODE Model', 'Gravity Model', 'Surface Tension Model', 'Mixed Model');
rlim([0, max(singular_radii_mm, [], 'omitnan') * 1.1]);
