%% plot_LD_polar.m
% Plot Aerodynamic Efficiency (L/D) vs Angle of Attack
clear; close all; clc;

% =====================================================================
% Global Plot Settings (LaTeX Interpreters)
% =====================================================================
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

% =====================================================================
% Load Data & Calculate
% =====================================================================
% Extract columns (assuming Col 1: Alpha, Col 2: CL, Col 3: CD)
polar_data = load('polar.txt');
alpha = polar_data(:, 1);
CL    = polar_data(:, 2);
CD    = polar_data(:, 3);

L_over_D = CL ./ CD;

% =====================================================================
% PLOT: L/D vs Alpha
% =====================================================================
figure('Name', 'L/D vs Alpha', 'Color', 'w', ...
       'Position', [100, 100, 500, 350]);
hold on; grid on;

% Plot the L/D curve (using the exact blue color from your reference script)
plot(alpha, L_over_D, '-', 'Color', [0.000, 0.447, 0.741], 'LineWidth', 1.5);

% Add the vertical line and text for Max L/D
% Note: Using LaTeX formatting for the label string
xline(4.5, '--', 'Max $L/D$ ($\alpha = 4.5^\circ$)', ...
    'Color', [0.850, 0.325, 0.098], ... % Using the red/orange from your reference
    'LineWidth', 1.2, ...
    'LabelHorizontalAlignment', 'right', ...
    'LabelVerticalAlignment', 'bottom', ...
    'FontSize', 10, ...
    'Interpreter', 'latex');

% Add Labels
xlabel('$\alpha \ [^\circ]$', 'FontSize', 13);
ylabel('$L/D \ [-]$', 'FontSize', 13);
title('Aero. Efficiency vs Angle of Attack (HQ-17 2D Airfoil)', 'FontSize', 14);
xlim([-4 12])
% Apply Axes Formatting (Box, LineWidth, Tick direction, and Axes FontSize)
set(gca, 'Box', 'on', 'LineWidth', 1, 'TickDir', 'in');
set(gca, 'FontSize', 11);