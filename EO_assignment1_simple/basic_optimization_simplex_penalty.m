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
    
    iter_count = 0;       
    total_aero_calls = 0; 
    
    %% ========================
    %  1. OPTIMIZATION SETTINGS (Back to 2D)
    %  ========================
    x_0 = [1.5, 0.75]; 
    lb  = [0.4, 0.2];
    ub  = [2.0, 0.8];
    max_iterations = 500; 
        
    %% ========================
    %  2. RUN CUSTOM SIMPLEX
    %  ========================
    fprintf('\n--- Starting Nelder-Mead Simplex (2D Penalty Method) ---\n');
    
    tol_f = 1e-4;
    tol_x = 1e-4;
    
    % Running your standard 2D simplex
    [c_root_opt, c_tip_opt] = simplex(@obj_fun, x_0, ub, lb, max_iterations, tol_f, tol_x);
    
    fprintf('\n======================================================\n');
    fprintf('OPTIMIZATION COMPLETED\n');
    fprintf('Total Simplex Iterations: %d\n', iter_count);
    fprintf('Total Aero Model Evaluations: %d\n', total_aero_calls);
    
    fprintf('Optimal Geometry : c_root = %5.4f | c_tip = %5.4f\n', c_root_opt, c_tip_opt);
    disp(history); 
    
    % --- SAVE TO .MAT FILE ---
    timestamp = datestr(now, 'yyyy-mm-dd_HHMM');
    filename = ['opt_penalty_history_p5', timestamp, '.mat']; % Renamed so it doesn't overwrite your nested run!
    
    save(filename, 'history', 'c_root_opt', 'c_tip_opt', 'total_aero_calls');
                   
    fprintf('History saved to: %s\n', filename);
    fprintf('======================================================\n\n');

    %% ========================
    %  3. SIMPLEX OBJECTIVE FUNCTION (2D PENALTY)
    %  ========================
    function f = obj_fun(x_bar)
        cr = x_bar(1);
        ct = x_bar(2);
        
        total_aero_calls = total_aero_calls + 1;
        
        % Constants
        V_fixed = 30.55; fixed_b2 = 7.5; fixed_twist = 0; h_cruise = 1500;
        
        % ---> THE FIXED ALPHA TRICK <---
        % We lock the aircraft at a realistic cruise angle (e.g., 4.5 degrees)
        alpha_fixed = 4.5; 
        
        % 1. Estimate Target Weight for this specific geometry
        aircraft = calc_planform(fixed_b2, cr, ct, fixed_twist);
        aero     = calc_atmos_properties(h_cruise, V_fixed, 'v', aircraft);
        [W_target, ~] = estimate_weight(aircraft, aero, V_fixed);
        
        % 2. Run Aero Model once 
        try
            [LD_val, L_tmp, ~] = run_model([fixed_b2, cr, ct, fixed_twist, V_fixed, alpha_fixed]);
            
            if isnan(L_tmp) 
                f = 1e6; 
                LD_val = NaN; err_pct = NaN;
            else 
                % 3. Calculate Lift Error
                err_pct = ((L_tmp - W_target) / W_target) * 100; 
                
                % 4. APPLY PENALTY
                % If Lift doesn't equal Weight, ruin the L/D score.
                % A multiplier of 5 means a 10% error adds 50 to the objective.
                penalty_multiplier = 5; 
                
                f = -LD_val + penalty_multiplier * abs(err_pct); 
            end
        catch
            f = 1e6;
            LD_val = NaN; err_pct = NaN;
        end
        
        % --- UPDATE HISTORY ---
        new_row = {iter_count, LD_val, cr, ct, err_pct, alpha_fixed};
        history = [history; new_row];
        
        fprintf('Iter %d: Root=%5.4f, Tip=%5.4f | L/D=%6.4f | Err=%+6.4f%% | Obj Fun=%6.4f\n', ...
                iter_count, cr, ct, LD_val, err_pct, f);
            
        iter_count = iter_count + 1; 
    end
end