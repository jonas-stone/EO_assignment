% script_4_sensitivity.m
% =========================================================================
% SENSITIVITY ANALYSIS — Quantitative Gradient Magnitudes for All 5 DVs
% =========================================================================
% PURPOSE: Quantify how much each variable influences the objective (L/D)
% and constraints. Complements the monotonicity analysis (qualitative:
% direction only) with quantitative gradient magnitudes.
%
% TWO METHODS IMPLEMENTED:
%   Option A: Average sensitivity from macro-sweep (script 2 data, recomputed)
%   Option B: Local central-difference gradient at the baseline point
%             using VALIDATED FD step sizes from noise analysis (Section 5.4)
%
% NORMALIZATION:
%   x_normalized = [1, 1, 1, 1, 0] at baseline
%   x_physical   = x_normalized .* x_baseline
%   x_baseline   = [1.5, 0.75, 0.4058, 30.55, 1]  (from constants())
%
% VARIABLE ORDERING:
%   Script:    [c_root, c_tip, alpha, V, twist]  → indices 1..5
%   run_model: [b2, c_root, c_tip, twist, V, alpha]
%   Mapping:   run_model([semi_span, x(1), x(2), x(5), x(4), x(3)])
% =========================================================================
clc; clear; close all;
addpath(genpath(pwd));

c_const = constants();
semi_span = 7.5;

% Normalized baseline
x_base = [1, 1, 1, 1, 0];
var_names = {'c_{root}', 'c_{tip}', '\alpha', 'V', 'twist'};
var_units = {'m', 'm', 'deg', 'm/s', 'deg'};

% Validated FD step sizes from noise analysis (Section 5.4 of porting doc)
% These are in NORMALIZED space (FinDiffRelStep values)
h_norm = [4.1e-3, 4.1e-3, 2.0e-2, 2.0e-3, 4.1e-2];

% Constraint parameters
CL_max      = 1.6;
V_stall_lim = 80 / 3.6;
W_limit     = 850 * c_const.g;

fprintf('=== SENSITIVITY ANALYSIS ===\n\n');

%% =====================================================================
%  OPTION B: LOCAL CENTRAL-DIFFERENCE GRADIENT AT BASELINE
%  (More rigorous — uses validated FD steps from noise analysis)
% ======================================================================
fprintf('--- Option B: Local Central-Difference Gradient at Baseline ---\n\n');

grad_LD     = zeros(1, 5);   % dLD/dx_norm
grad_stall  = zeros(1, 5);   % d(g_stall)/dx_norm
grad_weight = zeros(1, 5);   % d(g_weight)/dx_norm
grad_trim   = zeros(1, 5);   % d(h_trim)/dx_norm

for v = 1:5
    % Forward perturbation
    x_plus = x_base;
    x_plus(v) = x_base(v) + h_norm(v);
    x_phys_p = x_plus .* c_const.x_baseline;
    
    aircraft_p = calc_planform(semi_span, x_phys_p(1), x_phys_p(2), x_phys_p(5));
    aero_p = calc_atmos_properties(c_const.altitude, x_phys_p(4), 'v', aircraft_p);
    [W_p, ~] = estimate_weight(aircraft_p, aero_p, x_phys_p(4));
    [LD_p, L_p, ~] = run_model([semi_span, x_phys_p(1), x_phys_p(2), x_phys_p(5), x_phys_p(4), x_phys_p(3)]);
    
    CL_req_p = W_p / (0.5 * aero_p.rho * V_stall_lim^2 * aircraft_p.S);
    stall_p  = (CL_req_p / CL_max) - 1;
    wt_p     = (W_p / W_limit) - 1;
    trim_p   = (L_p - W_p) / W_p;
    
    % Backward perturbation
    x_minus = x_base;
    x_minus(v) = x_base(v) - h_norm(v);
    x_phys_m = x_minus .* c_const.x_baseline;
    
    aircraft_m = calc_planform(semi_span, x_phys_m(1), x_phys_m(2), x_phys_m(5));
    aero_m = calc_atmos_properties(c_const.altitude, x_phys_m(4), 'v', aircraft_m);
    [W_m, ~] = estimate_weight(aircraft_m, aero_m, x_phys_m(4));
    [LD_m, L_m, ~] = run_model([semi_span, x_phys_m(1), x_phys_m(2), x_phys_m(5), x_phys_m(4), x_phys_m(3)]);
    
    CL_req_m = W_m / (0.5 * aero_m.rho * V_stall_lim^2 * aircraft_m.S);
    stall_m  = (CL_req_m / CL_max) - 1;
    wt_m     = (W_m / W_limit) - 1;
    trim_m   = (L_m - W_m) / W_m;
    
    % Central difference
    grad_LD(v)     = (LD_p - LD_m) / (2 * h_norm(v));
    grad_stall(v)  = (stall_p - stall_m) / (2 * h_norm(v));
    grad_weight(v) = (wt_p - wt_m) / (2 * h_norm(v));
    grad_trim(v)   = (trim_p - trim_m) / (2 * h_norm(v));
    
    fprintf('  %s: h_norm = %.1e → LD+ = %.3f, LD- = %.3f → dLD/dx = %+.3f\n', ...
            var_names{v}, h_norm(v), LD_p, LD_m, grad_LD(v));
end

% Normalize gradients by the maximum for ranking
abs_grad_LD = abs(grad_LD);
max_grad = max(abs_grad_LD);
norm_sensitivity = abs_grad_LD / max_grad;

fprintf('\n--- Sensitivity Ranking (L/D, normalized) ---\n');
[sorted_sens, sort_idx] = sort(norm_sensitivity, 'descend');
fprintf('%-10s | %-12s | %-15s | %-15s\n', 'Rank', 'Variable', 'dLD/dx_norm', 'Rel. Sensitivity');
fprintf('%s\n', repmat('-', 1, 60));
for k = 1:5
    v = sort_idx(k);
    fprintf('  %d       | %-12s | %+12.4f    | %8.4f\n', ...
            k, var_names{v}, grad_LD(v), sorted_sens(k));
end

%% =====================================================================
%  OPTION A: AVERAGE SENSITIVITY FROM MACRO-SWEEP
%  (Coarser — uses full-range sweep like script 2)
% ======================================================================
fprintf('\n--- Option A: Average Macro-Sweep Sensitivity ---\n\n');

% Same physical ranges as script 2
ranges_phys = {
    linspace(0.6, 2.0, 15),   % c_root physical
    linspace(0.2, 0.8, 15),   % c_tip  physical
    linspace(-5,  8,   15),   % alpha  physical
    linspace(20,  55,  15),   % V      physical
    linspace(-6,  6,   15)    % twist  physical
};

avg_sens_LD = zeros(1, 5);

for v = 1:5
    sweep_phys = ranges_phys{v};
    n_pts = length(sweep_phys);
    LD_vec = NaN(1, n_pts);
    
    for i = 1:n_pts
        x_test = x_base;
        % Convert physical sweep value to normalized
        x_test(v) = sweep_phys(i) / c_const.x_baseline(v);
        x_phys = x_test .* c_const.x_baseline;
        
        try
            [LD_raw, ~, ~] = run_model( ...
                [semi_span, x_phys(1), x_phys(2), x_phys(5), x_phys(4), x_phys(3)]);
            LD_vec(i) = LD_raw;
        catch
            % leave NaN
        end
    end
    
    % Average sensitivity over the full range
    valid = ~isnan(LD_vec);
    if sum(valid) >= 2
        valid_idx = find(valid);
        avg_sens_LD(v) = (LD_vec(valid_idx(end)) - LD_vec(valid_idx(1))) / ...
                         (sweep_phys(valid_idx(end)) - sweep_phys(valid_idx(1)));
    end
    
    fprintf('  %s: dLD/d(phys) ≈ %+.4f per %s\n', var_names{v}, avg_sens_LD(v), var_units{v});
end

%% ===================== CONSTRAINT SENSITIVITY TABLE =====================
fprintf('\n--- Constraint Gradient Summary (at baseline, normalized) ---\n');
fprintf('%-10s | %-12s | %-12s | %-12s | %-12s\n', ...
        'Variable', 'dLD/dx', 'dStall/dx', 'dWeight/dx', 'dTrim/dx');
fprintf('%s\n', repmat('-', 1, 65));
for v = 1:5
    fprintf('%-10s | %+10.4f  | %+10.6f  | %+10.6f  | %+10.6f\n', ...
            var_names{v}, grad_LD(v), grad_stall(v), grad_weight(v), grad_trim(v));
end

%% ===================== VISUALIZATION =====================

figure('Name', 'Sensitivity Analysis', 'Position', [50 50 1200 500]);

% --- Panel 1: Bar chart of normalized L/D sensitivity ---
subplot(1, 2, 1);
bar_colors = [0.2 0.4 0.8; 0.2 0.4 0.8; 0.8 0.6 0.2; 0.2 0.4 0.8; 0.8 0.6 0.2];
b = bar(norm_sensitivity);
set(gca, 'XTickLabel', var_names);
ylabel('Normalized Sensitivity |dLD/dx|');
title('L/D Sensitivity Ranking (Local Gradient)');
grid on;
% Color bars: blue for strong, orange for weak
for k = 1:5
    if norm_sensitivity(k) < 0.1
        % Weak sensitivity — highlight
    end
end
% Add value labels on bars
for k = 1:5
    text(k, norm_sensitivity(k) + 0.03, sprintf('%.3f', norm_sensitivity(k)), ...
         'HorizontalAlignment', 'center', 'FontSize', 9);
end

% --- Panel 2: Constraint sensitivity heatmap ---
subplot(1, 2, 2);
constr_matrix = [grad_stall; grad_weight; grad_trim]';
imagesc(abs(constr_matrix'));
colorbar;
set(gca, 'XTick', 1:5, 'XTickLabel', var_names);
set(gca, 'YTick', 1:3, 'YTickLabel', {'Stall', 'Weight', 'Trim'});
title('Constraint Sensitivity Magnitude');
xlabel('Design Variable');
ylabel('Constraint');

sgtitle('Sensitivity Analysis: Objective and Constraints');

% Save figure
saveas(gcf, 'sensitivity_analysis.png');

%% ===================== CROSS-CHECK WITH NOISE =====================
fprintf('\n=== CROSS-CHECK: Sensitivity vs Noise ===\n');
fprintf('Noise floor from Q3D: ±0.015 L/D\n');
fprintf('Signal-to-noise ratio at validated FD step:\n\n');
noise_floor = 0.015;
for v = 1:5
    signal = abs(grad_LD(v) * h_norm(v));  % expected ΔLD from one FD step
    snr = signal / noise_floor;
    fprintf('  %s: signal = %.4f, SNR = %.2f %s\n', ...
            var_names{v}, signal, snr, ...
            ternary(snr > 3, '(GOOD)', ternary(snr > 1, '(MARGINAL)', '(POOR — noise dominates)')));
end

fprintf('\n=== CONCLUSIONS ===\n');
fprintf('Variables with strong sensitivity and good SNR → worth optimizing.\n');
fprintf('Variables with weak sensitivity or poor SNR → candidates for fixing.\n');
fprintf('  → Alpha near baseline is near the L/D optimum (self-optimizing).\n');
fprintf('  → Twist has very weak influence near zero.\n');
fprintf('  → c_root, c_tip, V are the primary optimization levers.\n');

% Helper function
function result = ternary(cond, if_true, if_false)
    if cond
        result = if_true;
    else
        result = if_false;
    end
end
