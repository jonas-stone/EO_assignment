% =========================================================================
% GLIDER WING OPTIMIZATION POST-PROCESSING (EQUAL-SIZED WINDOWS)
% =========================================================================
clear; clc; close all;

% --- 1. LOAD DATA ---
% Make sure to update this filename to your latest generated .mat file!
filename = 'opt_penalty_history_2026-04-14_1943.mat';
load(filename);

iter   = history.Iteration;
L_D    = history.L_D;
c_root = history.RootChord;
c_tip  = history.TipChord;
alpha  = history.Alpha;

% --- 2. CONSOLE OUTPUT ---
fprintf('\n======================================================\n');
fprintf('POST-PROCESSING RESULTS\n');
fprintf('======================================================\n');
if exist('total_aero_calls', 'var')
    fprintf('Total Aero Model Evaluations : %d\n', total_aero_calls);
end
fprintf('Baseline Geometry (Iter 0)   : c_root = %5.4f m | c_tip = %5.4f m\n', c_root(1), c_tip(1));
fprintf('Baseline L/D Ratio          : %6.4f\n', L_D(1));
fprintf('Baseline Trim Angle of Attack         : %5.4f deg\n', alpha(1));
fprintf('Optimal Geometry (Final)     : c_root = %5.4f m | c_tip = %5.4f m\n', c_root(end), c_tip(end));
fprintf('Maximized L/D Ratio          : %6.4f\n', L_D(end));
fprintf('Trim Angle of Attack         : %5.4f deg\n', alpha(end));
fprintf('======================================================\n\n');

% --- 3. PLOT SETTINGS ---
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

% Standard Window Size [width, height] in pixels
W = 700; 
H = 500;

% Base screen coordinates to cascade the windows
X0 = 100; Y0 = 300; offset = 30;

% Colors
col_LD = [0.00, 0.45, 0.74]; 
col_rt = [0.85, 0.33, 0.10]; 
col_tp = [0.93, 0.69, 0.13]; 
col_al = [0.47, 0.67, 0.19]; 

% =========================================================================
% FIGURE 1: L/D RATIO
% =========================================================================
figure('Name', 'L/D Convergence', 'Color', 'w', 'Position', [X0, Y0, W, H]);
plot(iter, L_D, '-o', 'LineWidth', 1.5, 'MarkerSize', 5, 'Color', col_LD, 'MarkerFaceColor', col_LD);
grid on;
xlabel('Iteration ID', 'FontSize', 12);
ylabel('Lift-to-Drag Ratio ($L/D$)', 'FontSize', 12);
title('\textbf{Aerodynamic Efficiency Convergence}', 'FontSize', 14);
set(gca, 'FontSize', 11, 'GridAlpha', 0.3);

% =========================================================================
% FIGURE 2: DESIGN VARIABLES (CHORDS)
% =========================================================================
figure('Name', 'Chord Convergence', 'Color', 'w', 'Position', [X0+offset, Y0-offset, W, H]);
plot(iter, c_root, '-s', 'LineWidth', 1.5, 'MarkerSize', 5, 'Color', col_rt, 'MarkerFaceColor', col_rt);
hold on; grid on;
plot(iter, c_tip, '-^', 'LineWidth', 1.5, 'MarkerSize', 5, 'Color', col_tp, 'MarkerFaceColor', col_tp);
xlabel('Iteration ID', 'FontSize', 12);
ylabel('Chord Length ($m$)', 'FontSize', 12);
title('\textbf{Design Variable Convergence}', 'FontSize', 14);
legend({'$c_{root}$', '$c_{tip}$'}, 'Location', 'northeast', 'FontSize', 11);
set(gca, 'FontSize', 11, 'GridAlpha', 0.3);

% =========================================================================
% FIGURE 3: TRIM ANGLE OF ATTACK
% =========================================================================
figure('Name', 'Alpha Convergence', 'Color', 'w', 'Position', [X0+2*offset, Y0-2*offset, W, H]);
plot(iter, alpha, '-d', 'LineWidth', 1.5, 'MarkerSize', 5, 'Color', col_al, 'MarkerFaceColor', col_al);
grid on;
xlabel('Iteration ID', 'FontSize', 12);
ylabel('Trim Angle of Attack $\alpha$ ($^\circ$)', 'FontSize', 12);
title('\textbf{Trim Angle Convergence}', 'FontSize', 14);
set(gca, 'FontSize', 11, 'GridAlpha', 0.3);

% =========================================================================
% FIGURE 4: WING PLANFORM
% =========================================================================
figure('Name', 'Wing Planform', 'Color', 'w', 'Position', [X0+3*offset, Y0-3*offset, W, H]);

b_half = 7.5; 
X_base = [-0.25*c_root(1), -0.25*c_tip(1), 0.75*c_tip(1), 0.75*c_root(1), -0.25*c_root(1)];
Y_span = [0, b_half, b_half, 0, 0];
X_opt  = [-0.25*c_root(end), -0.25*c_tip(end), 0.75*c_tip(end), 0.75*c_root(end), -0.25*c_root(end)];

plot(Y_span, X_base, '--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);
hold on; grid on;
fill(Y_span, X_opt, col_LD, 'FaceAlpha', 0.2, 'EdgeColor', col_LD, 'LineWidth', 2);
plot([0, b_half], [0, 0], '-.', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off'); 

xlabel('Spanwise Coordinate $y$ ($m$)', 'FontSize', 12);
ylabel('Chordwise Coordinate $x$ ($m$)', 'FontSize', 12);
title('\textbf{Half-Span Wing Planform Comparison}', 'FontSize', 14);
legend({'Baseline Wing', 'Optimized Wing'}, 'Location', 'southeast', 'FontSize', 11);

% --- THE FIX: Real aspect ratio and explicit boundaries ---
axis equal; % Ensures 1:1 physical ratio
set(gca, 'YDir', 'reverse', 'FontSize', 11, 'GridAlpha', 0.3);

% Force the horizontal axis to show the entire span with a 0.5m margin
xlim([-0.5, b_half + 0.5]); 

% Dynamically set the vertical limits based on the actual chord sizes 
max_c = max([c_root(1), c_root(end)]);
ylim([-0.5*max_c, 1.25*max_c]);

% =========================================================================
% FIGURE 5: LIFT DISTRIBUTION (BASELINE VS OPTIMIZED VS ELLIPTICAL)
% =========================================================================
figure('Name', 'Lift Distribution', 'Color', 'w', 'Position', [X0+4*offset, Y0-4*offset, W, H]);

% Ensure the data exists in the loaded .mat file
if exist('Y_opt', 'var') && exist('L_opt', 'var') && exist('W_target', 'var')
    
    % 1. Theoretical Elliptical Lift Distribution
    L_0 = (2 * W_target) / (pi * b_half);
    Y_ellipse = linspace(0, b_half, 100);
    L_ellipse = L_0 .* sqrt(1 - (Y_ellipse ./ b_half).^2);
    
    % Plot Ellipse (Dashed Black)
    plot(Y_ellipse, L_ellipse, '--k', 'LineWidth', 2);
    hold on; grid on;
    
    % 2. Plot Baseline Wing (Solid Gray with square markers)
    plot(Y_base, L_base, '-s', 'LineWidth', 1.5, 'MarkerSize', 5, ...
         'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5]);
         
    % 3. Plot Optimized Wing (Solid Blue with circle markers)
    plot(Y_opt, L_opt, '-o', 'LineWidth', 2, 'MarkerSize', 5, ...
         'Color', col_LD, 'MarkerFaceColor', col_LD);
    
    % Formatting
    xlabel('Spanwise Coordinate $y$ ($m$)', 'FontSize', 12);
    ylabel('Lift per unit span $L''(y)$ ($N/m$)', 'FontSize', 12);
    title('\textbf{Spanwise Lift Distribution Comparison}', 'FontSize', 14);
    
    legend({'Ideal Elliptical', 'Baseline Wing', 'Optimized Wing'}, ...
            'Location', 'southwest', 'FontSize', 11);
    
    set(gca, 'FontSize', 11, 'GridAlpha', 0.3);
    xlim([0, b_half + 0.2]);
    
    % Set dynamic Y-limits based on whichever curve is highest
    max_L = max([max(L_ellipse), max(L_base), max(L_opt)]);
    ylim([0, max_L * 1.15]); % 15% buffer at the top
else
    fprintf('Spanwise data not found in .mat file. Please re-run the updated optimizer.\n');
end