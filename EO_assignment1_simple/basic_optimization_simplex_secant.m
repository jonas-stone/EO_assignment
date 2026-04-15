function main()
    close all; clc;
    
    % --- Setup Paths ---
    projectRoot = pwd; 
    addpath(genpath(projectRoot));
    addpath(genpath("module_geometry_preprocess"));
    addpath(genpath("Q3D"));
    
    % --- Initialize History Storage ---
    history = table([], [], [], [], [], [], 'VariableNames', ...
              {'Iteration', 'L_D', 'RootChord', 'TipChord', 'Error_Pct', 'Alpha'});
    
    iter_count = 0;       % Counts Simplex steps
    total_aero_calls = 0; % Counts every single call to run_model (Secant + Final)
    
    % ---> MEMORY VARIABLES FOR THE EXTRA EVALUATION SHORTCUT <---
    last_eval_alpha = NaN;
    last_eval_LD    = NaN;
    %% ========================
    %  1. OPTIMIZATION SETTINGS
    %  ========================
    x_0 = [1.5, 0.75]; 
    lb  = [0.4, 0.2];
    ub  = [2.0, 0.8];
    max_iterations = 500; 
        
    %% ========================
    %  2. RUN CUSTOM SIMPLEX
    %  ========================
    fprintf('\n--- Starting Nelder-Mead Simplex Optimization ---\n');
    
    tol_f = 1e-2;
    tol_x = 1e-3;
    [c_root_opt, c_tip_opt] = simplex(@obj_fun, x_0, ub, lb, max_iterations, tol_f, tol_x);
    
    fprintf('\n======================================================\n');
    fprintf('OPTIMIZATION COMPLETED\n');
    fprintf('Total Simplex Iterations: %d\n', iter_count);
    fprintf('Total Aero Model Evaluations: %d\n', total_aero_calls);
    
    % Adjusted to 4 decimal places
    fprintf('Optimal Geometry : c_root = %5.4f | c_tip = %5.4f\n', c_root_opt, c_tip_opt);
    disp(history); 
    
    % =========================================================================
    % EXTRACT SPANWISE LIFT DISTRIBUTIONS (BASELINE VS OPTIMIZED)
    % =========================================================================
    fprintf('Extracting spanwise data for Baseline and Optimized wings...\n');
    
    % Constants used in your obj_fun
    V_fixed = 30.55; fixed_b2 = 7.5; fixed_twist = 0; h_cruise = 1500;
    
    % Setup dummy aircraft and aero structs to get local atmospheric properties
    aircraft_dummy = calc_planform(fixed_b2, c_root_opt, c_tip_opt, fixed_twist);
    aero_dummy     = calc_atmos_properties(h_cruise, V_fixed, 'v', aircraft_dummy);
    
    % Get dynamic pressure for dimensional lift (N/m) using local aero struct
    q_inf = 0.5 * aero_dummy.rho * V_fixed^2;
    
    % --- A. BASELINE WING (Trimmed at Iteration 0) ---
    c_root_base = x_0(1);
    c_tip_base  = x_0(2);
    alpha_base  = history.Alpha(1); % Grabs the trimmed alpha found on step 1
    
    [~, ~, base_Res] = run_model([fixed_b2, c_root_base, c_tip_base, fixed_twist, V_fixed, alpha_base]);
    Y_base   = base_Res.Section.Y(:)'; 
    Cl_base  = base_Res.Section.Cl(:)';
    cy_base  = interp1([0, fixed_b2], [c_root_base, c_tip_base], Y_base);
    L_base   = q_inf .* cy_base .* Cl_base; % Dimensional Lift in N/m
    
    % --- B. OPTIMIZED WING (Trimmed at Final Iteration) ---
    alpha_opt = history.Alpha(end);
    
    [~, ~, opt_Res] = run_model([fixed_b2, c_root_opt, c_tip_opt, fixed_twist, V_fixed, alpha_opt]);
    Y_opt   = opt_Res.Section.Y(:)'; 
    Cl_opt  = opt_Res.Section.Cl(:)';
    cy_opt  = interp1([0, fixed_b2], [c_root_opt, c_tip_opt], Y_opt);
    L_opt   = q_inf .* cy_opt .* Cl_opt; % Dimensional Lift in N/m
    
    % --- C. TARGET WEIGHT (For Elliptical Math) ---
    % Reuse the dummy structs created above
    [W_target, ~]  = estimate_weight(aircraft_dummy, aero_dummy, V_fixed);
    
    % --- SAVE EVERYTHING TO .MAT FILE ---
    timestamp = datestr(now, 'yyyy-mm-dd_HHMM');
    filename = ['opt_history_', timestamp, '.mat'];
    
    save(filename, 'history', 'c_root_opt', 'c_tip_opt', 'total_aero_calls', ...
                   'Y_base', 'L_base', 'Y_opt', 'L_opt', 'W_target');
                   
    fprintf('History and spanwise data saved to: %s\n', filename);
    fprintf('======================================================\n\n');

    %% ========================
    %  3. INTERNAL PHYSICS ENGINE
    %  =======================
    function [best_LD, best_alpha, final_err_pct, success] = trim_engine(c_root, c_tip, b2, twist, V, h)
        original_dir = pwd; 
        cleanupObj = onCleanup(@() cd(original_dir));
        
        aircraft = calc_planform(b2, c_root, c_tip, twist);
        aero     = calc_atmos_properties(h, V, 'v', aircraft);
        [W_target, ~] = estimate_weight(aircraft, aero, V);
        
        a1 = 0; a2 = 3; 
        tol_secant = 1e-3; max_iter = 15; 
        
        % Secant will call lift_residual multiple times
        [final_res, best_alpha] = secant(@(a) lift_residual(a, b2, c_root, c_tip, twist, V, W_target), ...
                                        a1, a2, tol_secant, max_iter, 6, -8);
        
        final_err_pct = final_res * 100;
        
        if abs(final_res) < 0.01 
            % ---> USE THE SAVED EVALUATION SHORTCUT <---
            % If the angle Secant returned is the exact one we just evaluated...
            if abs(best_alpha - last_eval_alpha) < 1e-6
                best_LD = last_eval_LD; % ...just grab the saved L/D!
                success = ~isnan(best_LD);
            else
                % Fallback: Run model again if secant somehow returned an older alpha
                try
                    [best_LD, ~, ~] = run_model([b2, c_root, c_tip, twist, V, best_alpha]);
                    total_aero_calls = total_aero_calls + 1; 
                    success = ~isnan(best_LD);
                catch
                    success = false; best_LD = 0;
                end
            end
        else
            success = false; best_LD = 0;
        end
    end
    % --- RAW RESIDUAL FUNCTION ---
    function res = lift_residual(a, b2, cr, ct, tw, V, W_t)
        try
            total_aero_calls = total_aero_calls + 1;
            
            % ---> CAPTURE BOTH L/D (LD_tmp) AND LIFT (L_tmp) <---
            [LD_tmp, L_tmp, ~] = run_model([b2, cr, ct, tw, V, a]);
            
            if isnan(L_tmp) 
                res = 1e6; 
            else 
                res = (L_tmp - W_t)/W_t; 
                
                % ---> UPDATE THE SHARED MEMORY <---
                last_eval_alpha = a;
                last_eval_LD    = LD_tmp;
            end
        catch
            res = 1e6;
        end
    end
    %% ========================
    %  4. SIMPLEX OBJECTIVE FUNCTION
    %  ========================
    function f = obj_fun(x_bar)
        cr = x_bar(1);
        ct = x_bar(2);
        
        V_fixed = 30.55; fixed_b2 = 7.5; fixed_twist = 0; h_cruise = 1500;
        
        [LD_val, alpha_val, err_pct, success] = trim_engine(cr, ct, fixed_b2, fixed_twist, V_fixed, h_cruise);
        
        if success
            f = -LD_val;
        else
            f = 1e6;
            LD_val = NaN; alpha_val = NaN; err_pct = NaN;
        end
        
        % --- UPDATE HISTORY ---
        new_row = {iter_count, LD_val, cr, ct, err_pct, alpha_val};
        history = [history; new_row];
        
        % ---> UPDATED CONSOLE OUTPUT TO 4 DECIMAL PLACES <---
        fprintf('Iter %d: Root=%5.4f, Tip=%5.4f | L/D=%6.4f | Alpha=%5.4f | Err=%+6.4f%% | AeroCalls=%d\n', ...
                iter_count, cr, ct, LD_val, alpha_val, err_pct, total_aero_calls);
            
        iter_count = iter_count + 1; 
    end
end