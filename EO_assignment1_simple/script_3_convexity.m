% script_3_convexity.m
% =========================================================================
% CONVEXITY CHECK — 2D Feasible Domain on the L=W Trim Surface
% =========================================================================
% PURPOSE: Map the objective (L/D) and stress constraint over the
% (c_root, c_tip) design space, with alpha SOLVED at each point so that
% the trim condition L = W is exactly satisfied.
%
% This gives the TRUE feasible-domain landscape that the optimizer will
% traverse, rather than an arbitrary fixed-alpha slice.
%
% WHAT TO LOOK FOR:
%   - Single peak with smooth contours → likely convex, one SQP run is fine
%   - Multiple peaks / ridges → non-convex, need multi-start or global search
%   - Stress boundary (red) cuts into a concave shape → non-convex feasible set
% =========================================================================
clc; clear; close all;
addpath(genpath(pwd));

c_const = constants();
semi_span = 7.5;

% Fixed conditions for the 2D slice
V_fixed     = 30.55;   % Cruise velocity [m/s]
twist_fixed = 0;        % Tip twist [deg]

% Grid resolution (15x15 = 225 points, each with ~5-10 fzero iterations)
n_root = 15;
n_tip  = 15;
root_vec = linspace(0.6, 2.0, n_root);
tip_vec  = linspace(0.2, 0.6, n_tip);
[Root_Grid, Tip_Grid] = meshgrid(root_vec, tip_vec);

% Pre-allocate output grids
LD_Grid     = NaN(size(Root_Grid));
Alpha_Grid  = NaN(size(Root_Grid));
Stress_Grid = NaN(size(Root_Grid));

sigma_allow = c_const.material.Fc_y / 1e6;  % [MPa]

% fzero options: loose tolerance is fine since we just need trim, not FD precision
fzero_opts = optimset('TolX', 0.01, 'Display', 'off');

% Initial guess for alpha trim — will be warm-started across the grid
alpha_guess = 2.0;

fprintf('=== 2D Trimmed Convexity Grid (%dx%d = %d points) ===\n', ...
        n_root, n_tip, n_root * n_tip);
fprintf('Each point solves fzero on alpha to satisfy L = W.\n');
fprintf('This will take a while...\n\n');

tic;
for i = 1:size(Root_Grid, 1)
    for j = 1:size(Root_Grid, 2)
        c_r = Root_Grid(i,j);
        c_t = Tip_Grid(i,j);
        
        % Compute target weight for this geometry
        aircraft = calc_planform(semi_span, c_r, c_t, twist_fixed);
        aero = calc_atmos_properties(c_const.altitude, V_fixed, 'v', aircraft);
        [W_target, ~] = estimate_weight(aircraft, aero, V_fixed);
        
        % Define trim residual: L(alpha) - W = 0
        trim_fun = @(a) get_lift_for_trim(semi_span, c_r, c_t, twist_fixed, V_fixed, a) - W_target;
        
        try
            % Solve for the alpha that trims the aircraft
            alpha_trim = fzero(trim_fun, alpha_guess, fzero_opts);
            
            % Evaluate full model (aero + structural) at trimmed alpha
            [LD_raw, L_check, stress_tmp] = run_model( ...
                [semi_span, c_r, c_t, twist_fixed, V_fixed, alpha_trim]);
            
            % Only accept if trim actually converged (2% tolerance)
            if ~isnan(LD_raw) && abs((L_check - W_target)/W_target) < 0.02
                LD_Grid(i,j)    = LD_raw;
                Alpha_Grid(i,j) = alpha_trim;
                
                if ~isnan(stress_tmp)
                    Stress_Grid(i,j) = (stress_tmp / sigma_allow) - 1;
                end
            end
            
            % Warm-start: use this solution as the guess for the next point
            alpha_guess = alpha_trim;
            
        catch
            % fzero failed or Q3D diverged — leave as NaN
        end
    end
    fprintf('  Row %d/%d complete.\n', i, size(Root_Grid, 1));
end
elapsed = toc;
fprintf('\nGrid complete in %.1f seconds.\n', elapsed);

%% ===================== PLOTTING =====================

figure('Name', 'Trimmed Convexity Map', 'Position', [50 50 1400 500]);

% --- Panel 1: L/D landscape with stress boundary ---
subplot(1, 3, 1);
contourf(Root_Grid, Tip_Grid, LD_Grid, 20, 'LineStyle', 'none');
cb1 = colorbar; ylabel(cb1, 'L/D');
hold on;
% Overlay stress constraint boundary (g = 0 line)
contour(Root_Grid, Tip_Grid, Stress_Grid, [0 0], 'r-', 'LineWidth', 3);
title('L/D on Trim Surface');
xlabel('Root Chord (m)'); ylabel('Tip Chord (m)');
legend('', 'Stress Limit (g=0)', 'Location', 'best');
grid on;

% --- Panel 2: Trim alpha required at each point ---
subplot(1, 3, 2);
contourf(Root_Grid, Tip_Grid, Alpha_Grid, 20, 'LineStyle', 'none');
cb2 = colorbar; ylabel(cb2, '\alpha_{trim} (deg)');
title('Trim Angle of Attack');
xlabel('Root Chord (m)'); ylabel('Tip Chord (m)');
grid on;

% --- Panel 3: Stress margin map (feasibility) ---
subplot(1, 3, 3);
contourf(Root_Grid, Tip_Grid, Stress_Grid, 20, 'LineStyle', 'none');
cb3 = colorbar; ylabel(cb3, 'g_{stress}  (\leq 0 feasible)');
hold on;
contour(Root_Grid, Tip_Grid, Stress_Grid, [0 0], 'r-', 'LineWidth', 3);
title('Stress Constraint Margin');
xlabel('Root Chord (m)'); ylabel('Tip Chord (m)');
grid on;

sgtitle('Convexity Check: Trimmed (L=W) Feasible Domain');

fprintf('\n=== Interpretation Guide ===\n');
fprintf('Panel 1: If L/D contours have a SINGLE peak → likely convex.\n');
fprintf('         Multiple peaks or saddle points → non-convex, multi-start needed.\n');
fprintf('Panel 2: Smooth alpha map → trim is well-behaved.\n');
fprintf('         Discontinuities → possible bifurcation in Q3D solver.\n');
fprintf('Panel 3: Red line = structural boundary. If it carves a concave\n');
fprintf('         indentation into the feasible region → non-convex domain.\n');

% =========================================================================
% LOCAL FUNCTION: Extract lift from run_model for fzero
% =========================================================================
function L = get_lift_for_trim(b2, c_root, c_tip, twist, V, alpha)
    [~, L, ~] = run_model([b2, c_root, c_tip, twist, V, alpha]);
    if isnan(L)
        L = 0;  % Return zero so fzero doesn't crash; the point will be
                 % filtered out by the 2% tolerance check above.
    end
end
