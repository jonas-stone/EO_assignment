function [x_opt, fval] = main()
    close all; clc;
    
    % Get the current folder where THIS script is saved
    projectRoot = pwd; 
    addpath(genpath(projectRoot));
    addpath(genpath("module_geometry_preprocess"));
    addpath(genpath("Q3D"));
    fprintf('Path reset. Project Root: %s\n', projectRoot);
    
    %% ========================
    %  GLOBAL CONSTANTS (Shared with nested functions)
    %  ========================
    % Define these exactly ONCE so physics always match
    c_const     = constants();
    h_cruise    = c_const.altitude;
    V_fixed     = 33.33;   % Kept at 27.78 for both objective and constraint!
    fixed_b2    = 7.5;     % 15m total span
    fixed_twist = 0;

    %% ========================
    %  INITIAL POINT & BOUNDS 
    %  ========================
    % Variables: [c_root, c_tip, alpha]
    x_0 = [1.5,  0.75,  2.0]; 
    lb  = [0.3,  0.2,  -4.0]; 
    ub  = [2.0,  0.8,   6.0]; 
        
    %% ========================
    %  OPTIMIZATION SETUP
    %  ========================
    options = optimoptions('fmincon', ...
        'Algorithm',            'sqp', ...
        'Display',              'iter', ...
        'FiniteDifferenceStepSize', 1e-2, ...
        'ConstraintTolerance',  5e-3, ...
        'OutputFcn',            @my_output_fun, ...
        'MaxIterations',        100, ...
        'MaxFunctionEvaluations', 1000, ...
        'StepTolerance',        1e-4, ...
        'OptimalityTolerance',  1e-5);
        
    %% ========================
    %  RUN FMINCON
    %  ========================
    fprintf('\n--- Starting SAND optimization ---\n');
    [x_opt, fval] = fmincon(@obj_fun, x_0, [], [], [], [], lb, ub, @nonlcon_fun, options);
    
    %% ========================
    %  FMINCON OBJECTIVE FUNCTION
    %  ========================
    function f = obj_fun(x_bar)
        c_root = x_bar(1);
        c_tip  = x_bar(2);
        alpha  = x_bar(3); 
    
        original_dir = pwd; 
        cleanupObj = onCleanup(@() cd(original_dir));
        
        try
            [LD_raw, ~, ~] = run_model([fixed_b2, c_root, c_tip, fixed_twist, V_fixed, alpha]);
            f = -LD_raw;
        catch
            % PENALTY: Do not use NaN. 1000 acts as an L/D of -1000 (Terrible)
            f = 1000; 
        end
    end

    %% ========================
    %  FMINCON NONLINEAR CONSTRAINT (L = W +/- 0.5%)
    %  ========================
    function [c, ceq] = nonlcon_fun(x_bar)
        % fmincon MUST receive two outputs!
        c = []; % No inequality constraints
        
        c_root = x_bar(1);
        c_tip  = x_bar(2);
        alpha  = x_bar(3);
        
        original_dir = pwd; 
        cleanupObj = onCleanup(@() cd(original_dir));
        
        % 1. Calculate Target Weight for this geometry
        aircraft = calc_planform(fixed_b2, c_root, c_tip, fixed_twist);
        aero     = calc_atmos_properties(h_cruise, V_fixed, 'v', aircraft);
        [W_target, ~] = estimate_weight(aircraft, aero, V_fixed);
        
        % 2. Calculate Actual Lift at the guessed alpha
        try
            [~, L_tmp, ~] = run_model([fixed_b2, c_root, c_tip, fixed_twist, V_fixed, alpha]);
            
            % Raw fractional error
            err = (L_tmp - W_target) / W_target; 
            ceq = err;
            
            fprintf('    Alpha: %5.2f | Raw L_tmp: %8.2f | Target W: %8.2f | Err: %+6.2f%%\n', ...
                    alpha, L_tmp, W_target, err*100);
        catch
            % PENALTY: Massive constraint violation if it crashes
            ceq = 1e6;
            fprintf('    Alpha: %5.2f | Raw L_tmp:  CRASHED | Target W: %8.2f\n', alpha, W_target);
        end
    end

    %% ========================
    %  CUSTOM OUTPUT FUNCTION
    %  ========================
    function stop = my_output_fun(x, optimValues, state)
        stop = false; 
        
        if isequal(state, 'iter')
            fprintf('\n======================================================\n');
            fprintf('ITERATION %d COMPLETED\n', optimValues.iteration);
            fprintf('Current Best L/D : %6.2f\n', -optimValues.fval); 
            fprintf('Current Geometry : c_root = %5.3f | c_tip = %5.3f | alpha = %5.2f deg\n', x(1), x(2), x(3));
            fprintf('Constraint Violation : %6.4f\n', optimValues.constrviolation);
            fprintf('======================================================\n\n');
        end
    end
end