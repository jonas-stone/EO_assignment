% script_3_convexity.m
% =========================================================================
% CONVEXITY CHECK — 2D Feasible Domain on the L=W Trim Surface
% =========================================================================
% PURPOSE: Map L/D over (c_root, c_tip) with alpha SOLVED via fzero at
% each point to satisfy L = W exactly. This shows the objective on the
% actual constraint surface the optimizer will traverse.
%
% NORMALIZATION:
%   x_normalized = [1, 1, 1, 1, 0] at baseline
%   x_physical   = x_normalized .* x_baseline
%   x_baseline   = [1.5, 0.75, 0.4058, 30.55, 1]  (from constants())
%
% CONSTRAINT BOUNDARIES OVERLAID (simplified model):
%   - Stall:  CL_req(V_stall) / CL_max = 1  (red)
%   - Weight: W_total / W_limit = 1          (magenta)
%
% WHAT TO LOOK FOR:
%   Single peak + smooth contours → likely convex, single SQP run sufficient.
%   Multiple peaks / saddle points → non-convex, need multi-start.
%   Constraint boundary carves concave shape → non-convex feasible domain.
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

% Fixed conditions for the 2D slice (normalized)
V_norm     = 1;       % normalized → physical = 30.55 m/s
twist_norm = 0;       % normalized → physical = 0 deg
V_fixed     = V_norm * c_const.x_baseline(4);      % 30.55 m/s
twist_fixed = twist_norm * c_const.x_baseline(5);   % 0 deg

% Grid in NORMALIZED space
% Physical range: c_root [0.6, 2.0] → normalized [0.4, 1.333]
% Physical range: c_tip  [0.2, 0.6] → normalized [0.267, 0.8]
n_root = 15;
n_tip  = 15;
root_norm_vec = linspace(0.6/c_const.x_baseline(1), 2.0/c_const.x_baseline(1), n_root);
tip_norm_vec  = linspace(0.2/c_const.x_baseline(2), 0.6/c_const.x_baseline(2), n_tip);

% Physical grid for plotting
root_phys_vec = root_norm_vec * c_const.x_baseline(1);
tip_phys_vec  = tip_norm_vec * c_const.x_baseline(2);
[Root_Grid, Tip_Grid] = meshgrid(root_phys_vec, tip_phys_vec);

% Normalized grid for internal use
[Root_Norm_Grid, Tip_Norm_Grid] = meshgrid(root_norm_vec, tip_norm_vec);

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
alpha_guess = 2.0;  % warm-started across grid (physical degrees)

fprintf('=== 2D Trimmed Convexity Grid (%dx%d = %d points) ===\n', ...
        n_root, n_tip, n_root * n_tip);
fprintf('Sweeping c_root and c_tip in normalized space.\n');
fprintf('Each point solves fzero on alpha (physical) for L = W.\n\n');

tic;
for i = 1:size(Root_Grid, 1)
    for j = 1:size(Root_Grid, 2)
        % Physical values from normalized grid
        c_r = Root_Grid(i,j);     % already physical (denormalized above)
        c_t = Tip_Grid(i,j);      % already physical
        
        % Target weight for this geometry
        aircraft = calc_planform(semi_span, c_r, c_t, twist_fixed);
        aero = calc_atmos_properties(c_const.altitude, V_fixed, 'v', aircraft);
        [W_target, ~] = estimate_weight(aircraft, aero, V_fixed);
        
        % Constraint margins that only depend on geometry (no alpha needed)
        CL_req = W_target / (0.5 * aero.rho * V_stall_lim^2 * aircraft.S);
        Stall_Grid(i,j) = (CL_req / CL_max) - 1;
        Wt_Grid(i,j)    = (W_target / W_limit) - 1;
        
        % Solve for trim alpha: L(alpha) = W
        % alpha here is in physical degrees (fzero operates in physical space)
        trim_fun = @(a) get_lift_for_trim(semi_span, c_r, c_t, twist_fixed, V_fixed, a) - W_target;
        
        try
            alpha_trim = fzero(trim_fun, alpha_guess, fzero_opts);
            
            % Evaluate at trimmed condition
            % run_model expects: [b2, c_root, c_tip, twist, V, alpha]
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

% Save figure
saveas(gcf, 'convexity_map.png');

%% ===================== CONVEXITY DIAGNOSTIC =====================
fprintf('\n=== Convexity Diagnostic ===\n');

% Count number of local maxima in the L/D grid
valid_LD = LD_Grid(~isnan(LD_Grid));
[max_LD, max_idx] = max(valid_LD);
fprintf('Global L/D maximum in grid: %.2f\n', max_LD);
fprintf('Number of valid grid points: %d / %d\n', sum(~isnan(LD_Grid(:))), numel(LD_Grid));

% Check for multiple peaks by looking at local maxima
% (a point higher than all 4 neighbors)
n_peaks = 0;
for i = 2:size(LD_Grid,1)-1
    for j = 2:size(LD_Grid,2)-1
        val = LD_Grid(i,j);
        if isnan(val); continue; end
        neighbors = [LD_Grid(i-1,j), LD_Grid(i+1,j), LD_Grid(i,j-1), LD_Grid(i,j+1)];
        if all(~isnan(neighbors)) && all(val > neighbors)
            n_peaks = n_peaks + 1;
            fprintf('  Local peak at c_root=%.2f m, c_tip=%.2f m: L/D=%.2f\n', ...
                    Root_Grid(i,j), Tip_Grid(i,j), val);
        end
    end
end

if n_peaks <= 1
    fprintf('ASSESSMENT: Single peak detected → LIKELY CONVEX.\n');
    fprintf('  Single SQP run from reasonable starting point should suffice.\n');
else
    fprintf('ASSESSMENT: %d peaks detected → POSSIBLY NON-CONVEX.\n', n_peaks);
    fprintf('  Multi-start strategy recommended.\n');
end

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
    % run_model expects: [b2, c_root, c_tip, twist, V, alpha]
    [~, L, ~] = run_model([b2, c_root, c_tip, twist, V, alpha]);
    if isnan(L)
        L = 0; % Prevents fzero crash; filtered by 2% check above
    end
end
