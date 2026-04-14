function main()
    close all; clc;
    
    % --- Setup Paths ---
    projectRoot = pwd; 
    addpath(genpath(projectRoot));
    addpath(genpath("module_geometry_preprocess"));
    addpath(genpath("Q3D"));
    
    fprintf('Path reset. Project Root: %s\n', projectRoot);

    %% ========================
    %  1. OPTIMIZATION SETTINGS
    %  ========================
    % Variables: [c_root, c_tip]
    x_0 = [1.5, 0.75]; 
    lb  = [0.4, 0.2];
    ub  = [2.0, 0.8];
    max_iterations = 60; 
        
    %% ========================
    %  2. RUN CUSTOM SIMPLEX
    %  ========================
    fprintf('\n--- Starting Nelder-Mead Simplex Optimization ---\n');
    
    % Calling your simplex.m file
    tol_simplex = 1e-5;
    [c_root_opt, c_tip_opt] = simplex(@obj_fun, x_0, ub, lb, max_iterations, tol_simplex);

    fprintf('\n======================================================\n');
    fprintf('OPTIMIZATION COMPLETED\n');
    fprintf('Optimal Geometry : c_root = %5.3f | c_tip = %5.3f\n', c_root_opt, c_tip_opt);
    
    % Final check of the result
    final_score = obj_fun([c_root_opt, c_tip_opt]);
    fprintf('Final Best L/D   : %6.2f\n', -final_score);
    fprintf('======================================================\n\n');

%% ========================
    %  3. INTERNAL PHYSICS ENGINE
    %  ========================
    function [best_LD, best_alpha, success] = trim_engine(c_root, c_tip, b2, twist, V, h)
        % Directory protection
        original_dir = pwd; 
        cleanupObj = onCleanup(@() cd(original_dir));

        aircraft = calc_planform(b2, c_root, c_tip, twist);
        aero     = calc_atmos_properties(h, V, 'v', aircraft);
        [W_target, ~] = estimate_weight(aircraft, aero, V);
        
        % --- SECANT METHOD ------ %
        a1 = 0; 
        a2 = 3; 
        tol_secant = 1e-3; 
        max_iter = 15; 
        max_limit = 6;
        low_limit = -8;
        
        % Call your custom root-finder
        [final_res, best_alpha] = secant(@(a) lift_residual(a, b2, c_root, c_tip, twist, V, W_target), a1, a2, tol_secant, max_iter, max_limit, low_limit);
        
        % If the final residual is essentially zero, we found trim!
        if abs(final_res) < 0.01 
            try
                [best_LD, ~, ~] = run_model([b2, c_root, c_tip, twist, V, best_alpha]);
                success = ~isnan(best_LD);
            catch
                success = false; best_LD = 0;
            end
        else
            success = false; best_LD = 0;
        end
    end

    % --- RAW RESIDUAL FUNCTION ---
    function res = lift_residual(a, b2, cr, ct, tw, V, W_t)
         % ---> ADD THIS TRACKER <---
        % This prints the angle the Secant method is about to test.
        fprintf('      [Secant] Testing Alpha = %6.3f deg... ', a);
        

        original_dir = pwd; 
        cleanupObj = onCleanup(@() cd(original_dir));
        try
            [~, L_tmp, ~] = run_model([b2, cr, ct, tw, V, a]);
            
            if isnan(L_tmp)
                res = 1e6; % Massive penalty for a stall/crash
                fprintf('FAILED (NaN)\n');
            else
                % RAW ERROR (Not squared!) 
                % If Lift > Weight, res is positive. If Lift < Weight, res is negative.
                res = (L_tmp - W_t)/W_t; 
                diff_pct = res * 100;
                fprintf('Lift = %7.2f N | Target = %7.2f N | Diff = %+6.2f%%\n', L_tmp, W_t, diff_pct);
            end
        catch
            res = 1e6;
            fprintf('CRASHED\n');
        end
    end

    %% ========================
    %  4. SIMPLEX OBJECTIVE FUNCTION
    %  ========================
    function f = obj_fun(x_bar)
        cr = x_bar(1);
        ct = x_bar(2);
        
        % Glider Flight Constants
        V_fixed     = 30.55; 
        fixed_b2    = 7.5;     
        fixed_twist = 0;
        h_cruise    = 1500; % Matches your previous constants

        % Call the trim engine
        [LD_val, ~, success] = trim_engine(cr, ct, fixed_b2, fixed_twist, V_fixed, h_cruise);
        
        if success
            f = -LD_val; % Minimize negative L/D
            fprintf('Eval: Root=%5.4f, Tip=%5.4f | L/D=%6.4f\n', cr, ct, LD_val);
        else
            f = 1e6; % Penalty for crashes
            fprintf('Eval: Root=%5.3f, Tip=%5.3f | FAILED\n', cr, ct);
        end
    end
end