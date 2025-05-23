% Script 1: Normalized h(theta) with LaTeX labels and mm units

r0 = 1.5e-3;
theta = linspace(0, 2*pi, 1000);
n_values = [3, 5, 7];

% Create figure with white background
figure('Color', 'w');

% Create polar axes
pax = polaraxes;
hold(pax, 'on');

% Loop through n values
for i = 1:length(n_values)
    n = n_values(i);

    % Gamma normalization
    G1 = gamma(1 + 1/n);
    G2 = gamma(1 + 2/n);
    R = r0 * sqrt((pi / 4) * (G2 / G1^2));

    % Superellipse definition
    rcos = abs(cos(theta)).^n;
    rsin = abs(sin(theta)).^n;
    r = R * (rcos + rsin).^(-1/n);

    dr_dtheta = gradient(r, theta);
    h = 0.5 * r.^2 ./ sqrt(r.^2 + dr_dtheta.^2);
    h_mm = h * 1000;  % Convert to mm
    h_norm = h_mm / max(h_mm);  % Normalize

    % Plot
    polarplot(pax, theta, h_norm, 'LineWidth', 1.5, ...
        'DisplayName', ['$n = ', num2str(n), '$']);
end

% Set title and legend
title(pax, 'Normalized $h(\theta)$ for Different $n$', 'Interpreter', 'latex');
legend(pax, 'Location', 'southoutside', 'Interpreter', 'latex');

% Set font and style
set(pax, 'FontName', 'Times', 'FontSize', 12);

% Format radial ticks to show mm units, including 0
% Get current radial tick labels
rticks = get(pax, 'RTickLabel');

% Change only the 0 to '0 mm'
rticks{1} = '0 mm';
set(pax, 'RTickLabel', rticks);


% Optional: increase line visibility
pax.LineWidth = 1;
