% script_2_monotonicity.m
% =========================================================================
% MONOTONICITY & CONSTRAINT BOUNDING — All 5 Design Variables
% =========================================================================
% PURPOSE: For each design variable, determine:
%   1. Does increasing it increase or decrease L/D? (monotonicity direction)
%   2. Which constraints bound it from going to +/-infinity?
%
% INTERPRETATION (Papalambros First Monotonicity Principle):
%   - Variable INCREASES L/D → must be bounded ABOVE by a constraint
%   - Variable DECREASES L/D → must be bounded BELOW by a constraint
%   - If no constraint crosses zero in the bounding direction, the problem
%     is ILL-POSED: fmincon will push to the box bound.
%
% CONSTRAINTS PLOTTED (g <= 0 formulation, negative = feasible):
%   1. Stress:  (sigma_crit / sigma_allow) - 1
%   2. Stall:   (CL_req_at_Vstall / CL_max) - 1
%   3. Weight:  (W_total / W_limit) - 1
%   4. Trim:    (L - W) / W  [equality, shown for reference]
% =========================================================================
clc; clear; close all;
addpath(genpath(pwd));

c_const = constants();
semi_span = 7.5;

% Baseline design point: [c_root, c_tip, alpha, V, twist]
x_base = [1.5, 0.75, -1.125, 30.55, 0];
var_names = {'c_{root} (m)', 'c_{tip} (m)', '\alpha (deg)', 'V (m/s)', 'twist (deg)'};

% Sweep ranges — chosen to span the physically meaningful domain
ranges = {
    linspace(0.6, 2.0, 30),  % c_root
    linspace(0.2, 0.8, 30),  % c_tip
    linspace(-5, 8, 30),     % alpha
    linspace(20, 55, 30),    % V
    linspace(-6, 6, 30)      % twist
};

% Constraint limits
sigma_allow  = c_const.material.Fc_y / 1e6;   % Compressive allowable [MPa]
CL_max       = 1.6;                            % Max section CL (airfoil)
V_stall_lim  = 80 / 3.6;                       % CS-22 stall speed limit [m/s]
W_limit      = 850 * c_const.g;                % FAI Open Class [N]

fprintf('=== Monotonicity Sweeps (5 Variables) ===\n');

for v = 1:5
    sweep = ranges{v};
    n_pts = length(sweep);
    
    LD_vec     = NaN(1, n_pts);
    stress_vec = NaN(1, n_pts);
    stall_vec  = NaN(1, n_pts);
    weight_vec = NaN(1, n_pts);
    trim_vec   = NaN(1, n_pts);
    
    fprintf('Sweeping %s (%d points)...\n', var_names{v}, n_pts);
    
    for i = 1:n_pts
        x_test = x_base;
        x_test(v) = sweep(i);
        
        % Build geometry and atmosphere (consistent with run_model internals)
        aircraft = calc_planform(semi_span, x_test(1), x_test(2), x_test(5));
        aero = calc_atmos_properties(c_const.altitude, x_test(4), 'v', aircraft);
        [W_target, ~] = estimate_weight(aircraft, aero, x_test(4));
        
        try
            [LD_raw, L_tmp, stress_tmp] = run_model( ...
                [semi_span, x_test(1), x_test(2), x_test(5), x_test(4), x_test(3)]);
            
            LD_vec(i) = LD_raw;
            
            % 1. Stress constraint: g = (sigma/sigma_allow) - 1 <= 0
            if ~isnan(stress_tmp)
                stress_vec(i) = (stress_tmp / sigma_allow) - 1;
            end
            
            % 2. Stall constraint: CL required at V_stall must not exceed CL_max
            %    This depends on geometry (S) and weight, NOT on cruise speed.
            CL_req = W_target / (0.5 * aero.rho * V_stall_lim^2 * aircraft.S);
            stall_vec(i) = (CL_req / CL_max) - 1;
            
            % 3. Weight constraint: g = (W/W_limit) - 1 <= 0
            weight_vec(i) = (W_target / W_limit) - 1;
            
            % 4. Trim (equality reference): (L - W) / W
            if ~isnan(L_tmp)
                trim_vec(i) = (L_tmp - W_target) / W_target;
            end
            
        catch ME
            fprintf('  Point %d failed: %s\n', i, ME.message);
        end
    end
    
    % --- Plot ---
    figure('Name', sprintf('Monotonicity — %s', var_names{v}), ...
           'Position', [50+40*v, 50+40*v, 950, 550]);
    
    % Left axis: Objective
    yyaxis left
    plot(sweep, LD_vec, 'b-', 'LineWidth', 2.5);
    ylabel('L/D  (Objective — maximize)');
    
    % Right axis: All constraint margins
    yyaxis right
    hold on;
    plot(sweep, stress_vec, 'r--',  'LineWidth', 1.5, 'DisplayName', 'Stress');
    plot(sweep, stall_vec,  'k-.',  'LineWidth', 1.5, 'DisplayName', 'Stall');
    plot(sweep, weight_vec, 'm:',   'LineWidth', 1.5, 'DisplayName', 'Weight');
    plot(sweep, trim_vec,   'c-',   'LineWidth', 1.0, 'DisplayName', 'Trim (L=W)');
    yline(0, 'g-', 'LineWidth', 2);
    ylabel('Constraint Margin  (g \leq 0 feasible)');
    hold off;
    
    title(sprintf('Monotonicity Sweep: %s', var_names{v}));
    xlabel(var_names{v});
    legend('L/D', 'Stress', 'Stall', 'Weight', 'Trim', 'Location', 'best');
    grid on;
end

fprintf('\n=== Interpretation Checklist ===\n');
fprintf('For EACH variable, answer:\n');
fprintf('  1. Is L/D monotonically increasing or decreasing?\n');
fprintf('  2. Does at least one constraint cross zero in the bounding direction?\n');
fprintf('     → YES: that constraint is CRITICAL (must be active at optimum).\n');
fprintf('     → NO:  the variable hits its box bound. Check if bounds are physical.\n');
fprintf('  3. Does the Trim curve cross zero? Where? That is the trimmed flight point.\n');
