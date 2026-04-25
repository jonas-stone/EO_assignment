function [x_opt, fval, eval_history, iter_history] = main()
    close all; clc;
    
    % --- Setup Paths ---
    projectRoot = pwd; 
    addpath(genpath(projectRoot));
    addpath(genpath("module_geometry_preprocess"));
    addpath(genpath("Q3D"));
    
    %% ========================
    %  1. SHARED DATA & HISTORY
    %  ========================
    % Constants
    c_const     = constants();
    h_cruise    = c_const.altitude;
    V_fixed     = 30.55;  
    fixed_b2    = 7.5;     
    fixed_twist = 0;

    % Global Counters & Tables
    eval_count = 0;
    eval_history = table([], [], [], [], [], [], 'VariableNames', ...
                   {'EvalNum', 'L_D', 'RootChord', 'TipChord', 'Error_Pct', 'Alpha'});
               
    iter_history = table([], [], [], [], [], [], 'VariableNames', ...
                   {'Iteration', 'L_D', 'RootChord', 'TipChord', 'Error_Pct', 'Alpha'});

    % Initialize Cache Memory
    cache.x    = NaN(1, 3); % [c_root, c_tip, alpha]
    cache.LD   = NaN;
    cache.L    = NaN;
    cache.W    = NaN;
    cache.err  = NaN;

    %% ========================
    %  2. OPTIMIZATION SETTINGS
    %  ========================
    x_0 = [0.8,  0.4,  0.0]; 
    lb  = [0.4,  0.2,  -4.0]; 
    ub  = [2.0,  0.8,   6.0]; 

    options = optimoptions('fmincon', ...
        'Algorithm',            'sqp', ...
        'Display',              'iter', ...
        'FiniteDifferenceStepSize', 1e-2, ...
        'ConstraintTolerance',  1e-3, ...
        'OutputFcn',            @my_output_fun, ...
        'MaxIterations',        100, ...
        'MaxFunctionEvaluations', 1000, ...
        'StepTolerance',        1e-4, ...
        'FunctionTolerance',    1e-4);
        
    %% ========================
    %  3. RUN FMINCON
    %  ========================
    fprintf('\n--- Starting Cached SQP Optimization ---\n');
    [x_opt, fval] = fmincon(@obj_fun, x_0, [], [], [], [], lb, ub, @nonlcon_fun, options);
    
    % Final Console Report
    fprintf('\n======================================================\n');
    fprintf('OPTIMIZATION COMPLETED\n');
    fprintf('Total Function Evaluations (Aero Calls): %d\n', eval_count);
    fprintf('Optimal Geometry : c_root = %5.4f | c_tip = %5.4f | alpha = %5.2f\n', x_opt(1), x_opt(2), x_opt(3));
    
    % Save data
    timestamp = datestr(now, 'yyyy-mm-dd_HHMM');
    save(['fmincon_results_', timestamp, '.mat'], 'eval_history', 'iter_history', 'x_opt');
    fprintf('Full history saved to .mat file.\n');

    %% ========================
    %  4. THE CACHE ENGINE
    %  ========================
    function update_cache(x)
        % Check if we already have this point in memory
        if norm(x - cache.x) < 1e-10
            return; 
        end
        
        cr = x(1); ct = x(2); alpha = x(3);
        
        % 1. Calculate Target Weight for this specific geometry
        aircraft = calc_planform(fixed_b2, cr, ct, fixed_twist);
        aero     = calc_atmos_properties(h_cruise, V_fixed, 'v', aircraft);
        [W_target, ~] = estimate_weight(aircraft, aero, V_fixed);
        
        original_dir = pwd; 
        cleanupObj = onCleanup(@() cd(original_dir));
        
        try
            % 2. Run Aero Model
            [LD_raw, L_tmp, ~] = run_model([fixed_b2, cr, ct, fixed_twist, V_fixed, alpha]);
            
            % 3. Update Cache
            cache.x   = x;
            cache.LD  = LD_raw;
            cache.L   = L_tmp;
            cache.W   = W_target;
            cache.err = (L_tmp - W_target) / W_target;
            
        catch
            % Failure Penalty
            cache.x   = x;
            cache.LD  = -1e6; 
            cache.L   = 0;
            cache.err = 1e6;
        end
        
        % 4. Record every single evaluation (Gradients + Steps)
        eval_count = eval_count + 1;
        new_eval = {eval_count, cache.LD, cr, ct, cache.err * 100, alpha};
        eval_history = [eval_history; new_eval];
    end

    %% ========================
    %  5. OBJ & CONSTRAINTS
    %  ========================
    function f = obj_fun(x)
        update_cache(x);
        f = -cache.LD;
    end

    function [c, ceq] = nonlcon_fun(x)
        update_cache(x);
        c   = []; 
        ceq = cache.err; % Equality constraint: L - W = 0
    end

    %% ========================
    %  6. OUTPUT FUNCTION (Iterations)
    %  ========================
    function stop = my_output_fun(x, optimValues, state)
        stop = false; 
        if isequal(state, 'iter')
            % When fmincon takes a step, record it
            new_iter = {optimValues.iteration, -optimValues.fval, x(1), x(2), cache.err * 100, x(3)};
            iter_history = [iter_history; new_iter];
            
            fprintf('  [ITER %d] L/D: %6.3f | Chord: %5.3f/%5.3f | Alpha: %5.2f | Err: %+6.3f%%\n', ...
                    optimValues.iteration, -optimValues.fval, x(1), x(2), x(3), cache.err*100);
        end
    end
end