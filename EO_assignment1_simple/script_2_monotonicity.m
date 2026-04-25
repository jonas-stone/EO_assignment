% script_2_monotonicity.m
% =========================================================================
% MONOTONICITY & CONSTRAINT BOUNDING — All 5 Design Variables
% =========================================================================
% PURPOSE: For each design variable, determine:
%   1. Does increasing it increase or decrease L/D?
%   2. Which active constraint bounds it from running to infinity?
%
% PAPALAMBROS FIRST MONOTONICITY PRINCIPLE:
%   Variable INCREASES L/D → must be bounded ABOVE by an active constraint.
%   Variable DECREASES L/D → must be bounded BELOW by an active constraint.
%   If no constraint crosses zero → variable hits box bound → check physics.
%
% ACTIVE CONSTRAINTS (simplified aero-only model):
%   1. Stall:   CL_req(V_stall) / CL_max  - 1  ≤ 0   (geometry + weight)
%   2. Weight:  W_total / W_limit  - 1          ≤ 0   (geometry only)
%   3. Trim:    (L - W) / W                     = 0   (equality, for reference)
% =========================================================================
clc; clear; close all;
addpath(genpath(pwd));

c_const = constants();
semi_span = 7.5;

% Baseline: [c_root, c_tip, alpha, V, twist]
x_base = [1.5, 0.75, -1.125, 30.55, 0];
var_names = {'c_{root} (m)', 'c_{tip} (m)', '\alpha (deg)', 'V (m/s)', 'twist (deg)'};

% Sweep ranges — physically meaningful domain for each variable
ranges = {
    linspace(0.6, 2.0, 30),   % c_root
    linspace(0.2, 0.8, 30),   % c_tip
    linspace(-5,  8,   30),   % alpha
    linspace(20,  55,  30),   % V
    linspace(-6,  6,   30)    % twist
};

% Constraint parameters
CL_max      = 1.6;             % Max section CL (airfoil limit)
V_stall_lim = 80 / 3.6;        % CS-22 stall speed [m/s]
W_limit     = 850 * c_const.g;  % FAI Open Class limit [N]

fprintf('=== Monotonicity Sweeps (5 Variables, Simplified Model) ===\n\n');

for v = 1:5
    sweep = ranges{v};
    n_pts = length(sweep);
    
    LD_vec    = NaN(1, n_pts);
    stall_vec = NaN(1, n_pts);
    wt_vec    = NaN(1, n_pts);
    trim_vec  = NaN(1, n_pts);
    
    fprintf('Sweeping %s (%d points)...\n', var_names{v}, n_pts);
    
    for i = 1:n_pts
        x_test = x_base;
        x_test(v) = sweep(i);
        
        aircraft = calc_planform(semi_span, x_test(1), x_test(2), x_test(5));
        aero = calc_atmos_properties(c_const.altitude, x_test(4), 'v', aircraft);
        [W_target, ~] = estimate_weight(aircraft, aero, x_test(4));
        
        try
            [LD_raw, L_tmp, ~] = run_model( ...
                [semi_span, x_test(1), x_test(2), x_test(5), x_test(4), x_test(3)]);
            
            LD_vec(i) = LD_raw;
            
            % Stall: CL required at V_stall must not exceed CL_max
            CL_req = W_target / (0.5 * aero.rho * V_stall_lim^2 * aircraft.S);
            stall_vec(i) = (CL_req / CL_max) - 1;
            
            % Weight
            wt_vec(i) = (W_target / W_limit) - 1;
            
            % Trim (equality reference)
            if ~isnan(L_tmp)
                trim_vec(i) = (L_tmp - W_target) / W_target;
            end
            
        catch ME
            fprintf('  Point %d failed: %s\n', i, ME.message);
        end
    end
    
    % --- Plot ---
    figure('Name', sprintf('Monotonicity — %s', var_names{v}), ...
           'Position', [50+40*v, 50+40*v, 900, 500]);
    
    yyaxis left
    plot(sweep, LD_vec, 'b-', 'LineWidth', 2.5);
    ylabel('L/D  (maximize)');
    
    yyaxis right
    hold on;
    plot(sweep, stall_vec, 'k-.',  'LineWidth', 1.5, 'DisplayName', 'Stall');
    plot(sweep, wt_vec,    'm:',   'LineWidth', 1.5, 'DisplayName', 'Weight');
    plot(sweep, trim_vec,  'c-',   'LineWidth', 1.0, 'DisplayName', 'Trim (L=W)');
    yline(0, 'g-', 'LineWidth', 2);
    ylabel('Constraint Margin  (g \leq 0 feasible)');
    hold off;
    
    title(sprintf('Monotonicity Sweep: %s', var_names{v}));
    xlabel(var_names{v});
    legend('L/D', 'Stall', 'Weight', 'Trim', 'Location', 'best');
    grid on;
end

fprintf('\n=== Interpretation Checklist ===\n');
fprintf('For EACH variable, answer:\n');
fprintf('  1. Is L/D monotonically increasing or decreasing over the range?\n');
fprintf('  2. Does at least one constraint cross g=0 in the bounding direction?\n');
fprintf('     YES → that constraint is CRITICAL (must be active at the optimum).\n');
fprintf('     NO  → variable hits its box bound. Are the bounds physical?\n');
fprintf('  3. Note: Stall and Weight depend only on geometry, not on alpha or V.\n');
fprintf('     They will appear flat on the alpha/V sweeps — that is correct.\n');
