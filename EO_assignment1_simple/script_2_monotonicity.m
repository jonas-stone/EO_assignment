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

% Normalized baseline: [c_root, c_tip, alpha, V, twist]
x_base = [1, 1, 1, 1, 0];
var_names = {'c_{root} (m)', 'c_{tip} (m)', '\alpha (deg)', 'V (m/s)', 'twist (deg)'};

% Sweep ranges in NORMALIZED space — converted from physical ranges
% Physical ranges:  c_root [0.6, 2.0], c_tip [0.2, 0.8], alpha [-5, 8],
%                   V [20, 55], twist [-6, 6]
% Normalized:       physical / x_baseline
ranges = {
    linspace(0.6/c_const.x_baseline(1), 2.0/c_const.x_baseline(1), 30),   % c_root norm
    linspace(0.2/c_const.x_baseline(2), 0.8/c_const.x_baseline(2), 30),   % c_tip  norm
    linspace(-5/c_const.x_baseline(3),  8/c_const.x_baseline(3),   30),   % alpha  norm
    linspace(20/c_const.x_baseline(4),  55/c_const.x_baseline(4),  30),   % V      norm
    linspace(-6/c_const.x_baseline(5),  6/c_const.x_baseline(5),   30)    % twist  norm
};

% Constraint parameters
CL_max      = 1.6;             % Max section CL (airfoil limit)
V_stall_lim = 80 / 3.6;        % CS-22 stall speed [m/s]
W_limit     = 850 * c_const.g;  % FAI Open Class limit [N]

fprintf('=== Monotonicity Sweeps (5 Variables, Normalized Space) ===\n\n');

% Store results for summary table
mono_direction = cell(1, 5);
bounding_constr = cell(1, 5);

for v = 1:5
    sweep = ranges{v};
    n_pts = length(sweep);
    
    LD_vec    = NaN(1, n_pts);
    stall_vec = NaN(1, n_pts);
    wt_vec    = NaN(1, n_pts);
    trim_vec  = NaN(1, n_pts);
    
    % Physical sweep values for axis labels
    sweep_phys = sweep * c_const.x_baseline(v);
    
    fprintf('Sweeping %s (%d points, normalized [%.2f, %.2f])...\n', ...
            var_names{v}, n_pts, sweep(1), sweep(end));
    
    for i = 1:n_pts
        x_test = x_base;
        x_test(v) = sweep(i);
        
        % Denormalize ALL variables to physical space
        x_phys = x_test .* c_const.x_baseline;
        
        aircraft = calc_planform(semi_span, x_phys(1), x_phys(2), x_phys(5));
        aero = calc_atmos_properties(c_const.altitude, x_phys(4), 'v', aircraft);
        [W_target, ~] = estimate_weight(aircraft, aero, x_phys(4));
        
        try
            [LD_raw, L_tmp, ~] = run_model( ...
                [semi_span, x_phys(1), x_phys(2), x_phys(5), x_phys(4), x_phys(3)]);
            
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
    
    % --- Determine monotonicity direction ---
    valid = ~isnan(LD_vec);
    if sum(valid) > 2
        p = polyfit(sweep_phys(valid), LD_vec(valid), 1);
        if p(1) > 0
            mono_direction{v} = 'INCREASING';
        elseif p(1) < 0
            mono_direction{v} = 'DECREASING';
        else
            mono_direction{v} = 'FLAT';
        end
    else
        mono_direction{v} = 'INSUFFICIENT DATA';
    end
    
    % --- Determine bounding constraint ---
    % Check which constraints cross g=0 within the sweep range
    stall_crosses = any(stall_vec(valid) < 0) && any(stall_vec(valid) > 0);
    wt_crosses    = any(wt_vec(valid) < 0) && any(wt_vec(valid) > 0);
    
    bounds = {};
    if stall_crosses; bounds{end+1} = 'Stall'; end
    if wt_crosses;    bounds{end+1} = 'Weight'; end
    if isempty(bounds)
        bounding_constr{v} = 'Box bound (no constraint crossing)';
    else
        bounding_constr{v} = strjoin(bounds, ' + ');
    end
    
    % --- Plot ---
    figure('Name', sprintf('Monotonicity — %s', var_names{v}), ...
           'Position', [50+40*v, 50+40*v, 900, 500]);
    
    yyaxis left
    plot(sweep_phys, LD_vec, 'b-', 'LineWidth', 2.5);
    ylabel('L/D  (maximize)');
    set(gca, 'YColor', 'b');
    
    yyaxis right
    hold on;
    plot(sweep_phys, stall_vec, 'r-.',  'LineWidth', 1.5, 'DisplayName', 'Stall (g_1)');
    plot(sweep_phys, wt_vec,    'm:',   'LineWidth', 1.5, 'DisplayName', 'Weight (g_2)');
    plot(sweep_phys, trim_vec,  'c-',   'LineWidth', 1.0, 'DisplayName', 'Trim (h_1)');
    yline(0, 'g-', 'LineWidth', 2, 'DisplayName', 'Feasibility boundary');
    ylabel('Constraint Margin  (g \leq 0 feasible)');
    set(gca, 'YColor', 'k');
    hold off;
    
    title(sprintf('Monotonicity: %s  [%s]', var_names{v}, mono_direction{v}));
    xlabel(var_names{v});
    legend('L/D', 'Stall (g_1)', 'Weight (g_2)', 'Trim (h_1)', 'g=0', 'Location', 'best');
    grid on;
    
    % Save figure
    saveas(gcf, sprintf('monotonicity_var%d.png', v));
    
    fprintf('  → Monotonicity: %s | Bounded by: %s\n\n', ...
            mono_direction{v}, bounding_constr{v});
end

%% ===================== SUMMARY TABLE =====================
fprintf('\n=== MONOTONICITY SUMMARY TABLE ===\n');
fprintf('%-16s | %-12s | %-35s | %-10s\n', ...
        'Variable', 'Direction', 'Bounding Constraint', 'Critical?');
fprintf('%s\n', repmat('-', 1, 80));
for v = 1:5
    % A constraint is critical if it must be active to bound the variable
    is_critical = ~contains(bounding_constr{v}, 'Box bound');
    fprintf('%-16s | %-12s | %-35s | %-10s\n', ...
            var_names{v}, mono_direction{v}, bounding_constr{v}, ...
            ternary(is_critical, 'YES', 'NO'));
end

fprintf('\n=== Interpretation Checklist ===\n');
fprintf('For EACH variable, answer:\n');
fprintf('  1. Is L/D monotonically increasing or decreasing over the range?\n');
fprintf('  2. Does at least one constraint cross g=0 in the bounding direction?\n');
fprintf('     YES → that constraint is CRITICAL (must be active at the optimum).\n');
fprintf('     NO  → variable hits its box bound. Are the bounds physical?\n');
fprintf('  3. Note: Stall and Weight depend only on geometry, not on alpha or V.\n');
fprintf('     They will appear flat on the alpha/V sweeps — that is correct.\n');
fprintf('  4. Alpha is bounded by trim (L=W): only one alpha trims the aircraft.\n');
fprintf('  5. V may be bounded only by box bounds — check if this is ill-posed.\n');

%% ===================== BOUNDEDNESS ANALYSIS =====================
fprintf('\n=== BOUNDEDNESS ANALYSIS ===\n');
fprintf('(Papalambros partial minimization — can f → -∞ for any variable?)\n\n');
for v = 1:5
    valid = ~isnan(LD_vec);
    fprintf('Variable %d (%s):\n', v, var_names{v});
    fprintf('  Monotonicity: %s\n', mono_direction{v});
    if strcmp(mono_direction{v}, 'INCREASING')
        fprintf('  L/D increases → needs UPPER bound.\n');
    elseif strcmp(mono_direction{v}, 'DECREASING')
        fprintf('  L/D decreases → needs LOWER bound.\n');
    else
        fprintf('  Self-bounded (interior optimum possible).\n');
    end
    fprintf('  Bounding mechanism: %s\n\n', bounding_constr{v});
end

% Helper function
function result = ternary(cond, if_true, if_false)
    if cond
        result = if_true;
    else
        result = if_false;
    end
end
