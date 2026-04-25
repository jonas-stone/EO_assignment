% script_1_noise.m
% =========================================================================
% NUMERICAL NOISE CHECK — All 5 Design Variables at fmincon FD Scale
% =========================================================================
% PURPOSE: Determine whether Q3D + structural model responses are smooth
% enough at the perturbation scale that fmincon's finite-difference
% gradient calculator will actually use (sqrt(eps)*max(1,|x|)).
%
% If the curves show flat plateaus, staircases, or random scatter,
% finite-difference gradients will be garbage and SQP will fail.
% =========================================================================
clc; clear; close all;
addpath(genpath(pwd));

c_const = constants();
semi_span = 7.5;

% Baseline design point
% Order: [c_root, c_tip, alpha, V, twist]
x_base = [1.5, 0.75, 0.4058, 30.55, 0];
var_names = {'c_{root} (m)', 'c_{tip} (m)', '\alpha (deg)', 'V (m/s)', 'twist (deg)'};

n_pts = 50;

% =========================================================================
% MAPPING NOTE:
% x_base index:    1=c_root  2=c_tip  3=alpha  4=V  5=twist
% run_model expects: [b2, c_root, c_tip, twist, V, alpha]
%                      1    2       3      4     5   6
% So the call is: run_model([semi_span, x(1), x(2), x(5), x(4), x(3)])
% =========================================================================

fprintf('=== Numerical Noise Check (fmincon FD scale) ===\n');

figure('Name', 'Noise Check — L/D', 'Position', [50 50 1400 800]);
figure('Name', 'Noise Check — (L-W)/W', 'Position', [100 100 1400 800]);
figure('Name', 'Noise Check — Stress', 'Position', [150 150 1400 800]);

for v = 1:5
    x_nominal = x_base(v);
    
    % This is the step size fmincon will use for central differences
    delta = sqrt(eps) * max(1, abs(x_nominal));
    
    sweep = linspace(x_nominal - delta, x_nominal + delta, n_pts);
    
    LD_vec  = zeros(1, n_pts);
    LW_err  = zeros(1, n_pts);
    stress_vec = zeros(1, n_pts);
    
    fprintf('Variable %d: %s | nominal = %.4f | delta = %.2e\n', ...
            v, var_names{v}, x_nominal, delta);
    
    for i = 1:n_pts
        x_test = x_base;
        x_test(v) = sweep(i);
        
        % Compute W_target externally (same function, same args as run_model)
        aircraft = calc_planform(semi_span, x_test(1), x_test(2), x_test(5));
        aero = calc_atmos_properties(c_const.altitude, x_test(4), 'v', aircraft);
        [W_target, ~] = estimate_weight(aircraft, aero, x_test(4));
        
        % Run the coupled model
        [LD_raw, L_tmp, stress_tmp] = run_model( ...
            [semi_span, x_test(1), x_test(2), x_test(5), x_test(4), x_test(3)]);
        
        LD_vec(i)     = LD_raw;
        stress_vec(i) = stress_tmp;
        
        if ~isnan(L_tmp) && W_target > 0
            LW_err(i) = (L_tmp - W_target) / W_target;
        else
            LW_err(i) = NaN;
        end
    end
    
    % --- Plot L/D ---
    figure(1);
    subplot(2, 3, v);
    plot(sweep, LD_vec, '-o', 'MarkerSize', 2, 'LineWidth', 1);
    xlabel(var_names{v}); ylabel('L/D'); grid on;
    title(sprintf('L/D vs %s  (\\Delta=%.1e)', var_names{v}, delta));
    
    % --- Plot (L-W)/W ---
    figure(2);
    subplot(2, 3, v);
    plot(sweep, LW_err, '-o', 'MarkerSize', 2, 'LineWidth', 1);
    xlabel(var_names{v}); ylabel('(L-W)/W'); grid on;
    title(sprintf('Trim Error vs %s', var_names{v}));
    
    % --- Plot Stress ---
    figure(3);
    subplot(2, 3, v);
    plot(sweep, stress_vec, '-o', 'MarkerSize', 2, 'LineWidth', 1);
    xlabel(var_names{v}); ylabel('\sigma_{crit} (MPa)'); grid on;
    title(sprintf('Stress vs %s', var_names{v}));
end

fprintf('\n=== Interpretation Guide ===\n');
fprintf('SMOOTH curve  → FD gradients will work for this variable.\n');
fprintf('STAIRCASE     → Q3D internal tolerance is coarser than FD step.\n');
fprintf('               Fix: increase FinDiffRelStep in fmincon options,\n');
fprintf('               or tighten AC.Aero.MaxIterIndex in Q3D.\n');
fprintf('RANDOM SCATTER→ Severe numerical noise. Consider derivative-free method.\n');
