% Script: q^(3/4)(theta) for different n, with experimental marker in legend

r0 = 1.5e-3;
theta = linspace(0, 2*pi, 1000);
U_0 = 2 / (6e4);
n_values = [3, 5, 7];

% Create polar plot with white background
figure('Color', 'w');
pax = polaraxes;
hold(pax, 'on');

for i = 1:length(n_values)
    n = n_values(i);

    % Gamma-based normalization
    G1 = gamma(1 + 1/n);
    G2 = gamma(1 + 2/n);
    R = r0 * sqrt((pi / 4) * (G2 / G1^2));

    % Superellipse radius
    rcos = abs(cos(theta)).^n;
    rsin = abs(sin(theta)).^n;
    r = R * (rcos + rsin).^(-1/n);

    % Derivatives and flow quantities
    dr_dtheta = gradient(r, theta);
    h = 0.5 * r.^2 ./ sqrt(r.^2 + dr_dtheta.^2);
    dl_dtheta = sqrt(r.^2 + dr_dtheta.^2);
    q = U_0 .* h .* dl_dtheta;
    q34 = 20e8 * q.^(3/4);  % Unnormalized scaled function

    % Plot model curve
    polarplot(pax, theta, q34, 'LineWidth', 1.5, ...
        'DisplayName', ['$n = ', num2str(n), '$']);
end

% Dummy marker for experimental data in legend
polarplot(pax, NaN, NaN, 'ko', 'MarkerFaceColor', 'k', ...
    'DisplayName', 'Experimental');

% Title and legend
title(pax, '$20 \times 10^8 \cdot q^{3/4}(\theta)$ with Experimental Data', ...
    'Interpreter', 'latex');
legend(pax, 'Location', 'southoutside', 'Interpreter', 'latex');
set(pax, 'FontName', 'Times', 'FontSize', 12);

% Optional: change radial tick for 0
rticks = get(pax, 'RTickLabel');
rticks{1} = '0 mm';
set(pax, 'RTickLabel', rticks);
