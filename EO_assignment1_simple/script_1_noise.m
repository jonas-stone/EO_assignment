% script_1_noise.m
% =========================================================================
% NUMERICAL NOISE CHECK — All 5 Design Variables at fmincon FD Scale
% =========================================================================
% PURPOSE: Determine whether Q3D responses are smooth at the perturbation
% scale that fmincon's finite-difference gradient will actually use.
%
% fmincon default FD step: sqrt(eps) * max(1, |x|)  ≈  1e-8 to 1e-7
%
% RESPONSES CHECKED (simplified aero-only model):
%   1. L/D  (objective)
%   2. (L - W) / W  (trim equality constraint)
%
% WHAT TO LOOK FOR:
%   SMOOTH monotonic line → FD gradients reliable for this variable.
%   FLAT PLATEAUS / STAIRCASES → solver tolerance coarser than FD step.
%       Fix: increase fmincon's FinDiffRelStep or tighten Q3D MaxIterIndex.
%   RANDOM SCATTER → severe noise. Consider derivative-free method.
% =========================================================================
clc; clear; close all;
addpath(genpath(pwd));

c_const = constants();
semi_span = 7.5;

% Baseline design point: [c_root, c_tip, alpha, V, twist]
x_base = [1, 1, 1, 1, 0];
var_names = {'c_{root} (m)', 'c_{tip} (m)', '\alpha (deg)', 'V (m/s)', 'twist (deg)'};

n_pts = 50;

% =========================================================================
% MAPPING:
%   x_base index:       1=c_root  2=c_tip  3=alpha  4=V     5=twist
%   run_model expects:  [b2,      c_root,  c_tip,   twist,  V,  alpha]
%   Call: run_model([semi_span, x(1), x(2), x(5), x(4), x(3)])
% =========================================================================

fprintf('=== Numerical Noise Check (fmincon FD scale) ===\n\n');

figure('Name', 'Noise Check — L/D',     'Position', [50  50  1400 700]);
figure('Name', 'Noise Check — Trim Error', 'Position', [100 100 1400 700]);
% Baseline design point: [c_root, c_tip, alpha, V, twist]
delta = [1e-1,1e-1,5e-1,5e-2,1];
for v = 1:5
    x_nominal = x_base(v);
    %delta = 5e-1;%sqrt(eps) * max(1, abs(x_nominal)); % fmincon's actual FD step
    sweep = linspace(x_nominal - delta(v), x_nominal + delta(v), n_pts);
    
    LD_vec = zeros(1, n_pts);
    LW_err = zeros(1, n_pts);
    
    fprintf('Variable %d: %-12s | nominal = %8.4f | delta = %.2e\n', ...
            v, var_names{v}, x_nominal, delta(v));
    
    for i = 1:n_pts
        x_test = x_base;
        x_test(v) = sweep(i);                        % perturb in normalized space
        x_phys = x_test.*c_const.x_baseline;       % denormalize ALL 5 variables
    
        aircraft = calc_planform(semi_span, x_phys(1), x_phys(2), x_phys(5));
        aero = calc_atmos_properties(c_const.altitude, x_phys(4), 'v', aircraft);
        [W_target, ~] = estimate_weight(aircraft, aero, x_phys(4));
    
        [LD_raw, L_tmp, ~] = run_model( ...
        [semi_span, x_phys(1), x_phys(2), x_phys(5), x_phys(4), x_phys(3)]);
        
        LD_vec(i) = LD_raw;
        if ~isnan(L_tmp) && W_target > 0
            LW_err(i) = (L_tmp - W_target) / W_target;
        else
            LW_err(i) = NaN;
        end
    end
    
    % --- L/D ---
    figure(1);
    subplot(2, 3, v);
    plot(sweep*c_const.x_baseline(v), LD_vec, '-o', 'MarkerSize', 2, 'LineWidth', 1);
    xlabel(var_names{v}); ylabel('L/D'); grid on;
    title(sprintf('L/D vs %s  (\\Delta = %.1e)', var_names{v}, delta(v)));
    
    % --- Trim error ---
    figure(2);
    subplot(2, 3, v);
    plot(sweep*c_const.x_baseline(v), LW_err, '-o', 'MarkerSize', 2, 'LineWidth', 1);
    xlabel(var_names{v}); ylabel('(L-W)/W'); grid on;
    title(sprintf('Trim Error vs %s', var_names{v}));
end

figure(1); sgtitle('Noise Check: Objective (L/D)');
figure(2); sgtitle('Noise Check: Equality Constraint (L-W)/W');
