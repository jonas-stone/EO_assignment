% script_3_convexity.m
% =========================================================================
% CONVEXITY CHECK — 2D Feasible Domain on the L=W Trim Surface
% =========================================================================
% PURPOSE: Map L/D over (c_root, c_tip) with alpha SOLVED via fzero at
% each point to satisfy L = W exactly. This shows the objective on the
% actual constraint surface the optimizer will traverse.
%
% CONSTRAINT BOUNDARIES OVERLAID (simplified model):
%   - Stall:  CL_req(V_stall) / CL_max = 1  (red)
%   - Weight: W_total / W_limit = 1          (magenta)
%
% WHAT TO LOOK FOR:
%   Single peak + smooth contours → likely convex, single SQP run sufficient.
%   Multiple peaks / saddle points → non-convex, need multi-start.
%   Constraint boundary carves concave shape → non-convex feasible domain.
% =========================================================================
clc; clear; close all;
addpath(genpath(pwd));

c_const = constants();
semi_span = 7.5;

% Fixed conditions for the 2D slice
V_fixed     = 30.55;
twist_fixed = 0;

% Grid
n_root = 15;
n_tip  = 15;
root_vec = linspace(0.6, 2.0, n_root);
tip_vec  = linspace(0.2, 0.6, n_tip);
[Root_Grid, Tip_Grid] = meshgrid(root_vec, tip_vec);

% Pre-allocate
LD_Grid    = NaN(size(Root_Grid));
Alpha_Grid = NaN(size(Root_Grid));
Stall_Grid = NaN(size(Root_Grid));
Wt_Grid    = NaN(size(Root_Grid));

% Constraint parameters
CL_max      = 1.6;
V_stall_lim = 80 / 3.6;
W_limit     = 850 * c_const.g;

fzero_opts = optimset('TolX', 0.01, 'Display', 'off');
alpha_guess = 2.0;  % warm-started across grid

fprintf('=== 2D Trimmed Convexity Grid (%dx%d = %d points) ===\n', ...
        n_root, n_tip, n_root * n_tip);
fprintf('Each point solves fzero on alpha for L = W.\n\n');

tic;
for i = 1:size(Root_Grid, 1)
    for j = 1:size(Root_Grid, 2)
        c_r = Root_Grid(i,j);
        c_t = Tip_Grid(i,j);
        
        % Target weight for this geometry
        aircraft = calc_planform(semi_span, c_r, c_t, twist_fixed);
        aero = calc_atmos_properties(c_const.altitude, V_fixed, 'v', aircraft);
        [W_target, ~] = estimate_weight(aircraft, aero, V_fixed);
        
        % Constraint margins that only depend on geometry (no alpha needed)
        CL_req = W_target / (0.5 * aero.rho * V_stall_lim^2 * aircraft.S);
        Stall_Grid(i,j) = (CL_req / CL_max) - 1;
        Wt_Grid(i,j)    = (W_target / W_limit) - 1;
        
        % Solve for trim alpha: L(alpha) = W
        trim_fun = @(a) get_lift_for_trim(semi_span, c_r, c_t, twist_fixed, V_fixed, a) - W_target;
        
        try
            alpha_trim = fzero(trim_fun, alpha_guess, fzero_opts);
            
            % Evaluate at trimmed condition
            [LD_raw, L_check, ~] = run_model( ...
                [semi_span, c_r, c_t, twist_fixed, V_fixed, alpha_trim]);
            
            % Accept only if trim actually converged
            if ~isnan(LD_raw) && abs((L_check - W_target)/W_target) < 0.02
                LD_Grid(i,j)    = LD_raw;
                Alpha_Grid(i,j) = alpha_trim;
            end
            
            alpha_guess = alpha_trim;  % warm-start
            
        catch
            % fzero or Q3D failed — leave NaN
        end
    end
    fprintf('  Row %d/%d complete.\n', i, size(Root_Grid, 1));
end
fprintf('\nGrid complete in %.1f seconds.\n', toc);

%% ===================== PLOTTING =====================

figure('Name', 'Trimmed Convexity Map', 'Position', [50 50 1200 500]);

% --- Panel 1: L/D landscape with constraint boundaries ---
subplot(1, 2, 1);
contourf(Root_Grid, Tip_Grid, LD_Grid, 20, 'LineStyle', 'none');
cb1 = colorbar; ylabel(cb1, 'L/D');
hold on;
contour(Root_Grid, Tip_Grid, Stall_Grid, [0 0], 'r-',  'LineWidth', 3);
contour(Root_Grid, Tip_Grid, Wt_Grid,    [0 0], 'm--', 'LineWidth', 3);
title('L/D on Trim Surface (L=W satisfied)');
xlabel('Root Chord (m)'); ylabel('Tip Chord (m)');
legend('', 'Stall Limit', 'Weight Limit', 'Location', 'best');
grid on;

% --- Panel 2: Trim alpha required ---
subplot(1, 2, 2);
contourf(Root_Grid, Tip_Grid, Alpha_Grid, 20, 'LineStyle', 'none');
cb2 = colorbar; ylabel(cb2, '\alpha_{trim} (deg)');
hold on;
contour(Root_Grid, Tip_Grid, Stall_Grid, [0 0], 'r-',  'LineWidth', 3);
contour(Root_Grid, Tip_Grid, Wt_Grid,    [0 0], 'm--', 'LineWidth', 3);
title('Trim Angle of Attack');
xlabel('Root Chord (m)'); ylabel('Tip Chord (m)');
legend('', 'Stall Limit', 'Weight Limit', 'Location', 'best');
grid on;

sgtitle('Convexity Check: Trimmed Feasible Domain');

fprintf('\n=== Interpretation ===\n');
fprintf('Panel 1: Single L/D peak inside feasible region → likely convex.\n');
fprintf('         Multiple peaks or saddle → non-convex, multi-start needed.\n');
fprintf('Panel 2: Smooth alpha field → trim is well-behaved.\n');
fprintf('         Jumps/discontinuities → possible Q3D bifurcation.\n');
fprintf('Feasible region = below/left of BOTH constraint lines.\n');

% =========================================================================
% LOCAL FUNCTION
% =========================================================================
function L = get_lift_for_trim(b2, c_root, c_tip, twist, V, alpha)
    [~, L, ~] = run_model([b2, c_root, c_tip, twist, V, alpha]);
    if isnan(L)
        L = 0; % Prevents fzero crash; filtered by 2% check above
    end
end
