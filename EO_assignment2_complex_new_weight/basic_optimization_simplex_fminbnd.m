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
    lb  = [0.8, 0.2];
    ub  = [2.0, 0.8];
    max_iterations = 60; 
        
    %% ========================
    %  2. RUN CUSTOM SIMPLEX
    %  ========================
    fprintf('\n--- Starting Nelder-Mead Simplex Optimization ---\n');
    
    % Calling your simplex.m file
    [c_root_opt, c_tip_opt] = simplex(@obj_fun, x_0, ub, lb, max_iterations);

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
        
        % Solve for Alpha (Trim)
        opts = optimset('Display', 'none', 'TolX', 1e-3);
        [best_alpha, min_err] = fminbnd(@(a) lift_err(a, b2, c_root, c_tip, twist, V, W_target), -4, 8, opts);
        
        if min_err < 0.001 
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

    function err = lift_err(a, b2, cr, ct, tw, V, W_t)
        original_dir = pwd; 
        cleanupObj = onCleanup(@() cd(original_dir));
        try
            [~, L_tmp, ~] = run_model([b2, cr, ct, tw, V, a]);
            if isnan(L_tmp), err = 1e6; else, err = ((L_tmp - W_t)/W_t)^2; end
        catch
            err = 1e6;
        end
    end

    %% ========================
    %  4. SIMPLEX OBJECTIVE FUNCTION
    %  ========================
    function f = obj_fun(x_bar)
        cr = x_bar(1);
        ct = x_bar(2);
        
        % Glider Flight Constants
        V_fixed     = 27.78; 
        fixed_b2    = 7.5;     
        fixed_twist = 0;
        h_cruise    = 1500; % Matches your previous constants

        % Call the trim engine
        [LD_val, ~, success] = trim_engine(cr, ct, fixed_b2, fixed_twist, V_fixed, h_cruise);
        
        if success
            f = -LD_val; % Minimize negative L/D
            fprintf('Eval: Root=%5.3f, Tip=%5.3f | L/D=%6.2f\n', cr, ct, LD_val);
        else
            f = 1e6; % Penalty for crashes
            fprintf('Eval: Root=%5.3f, Tip=%5.3f | FAILED\n', cr, ct);
        end
    end
end