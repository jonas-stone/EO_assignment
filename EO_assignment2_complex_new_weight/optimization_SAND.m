function [x_opt, fval] = main()
close all; clc;
    
    % Get the current folder where THIS script is saved
    projectRoot = pwd; 
    
    % Add EVERYTHING in the project to the path recursively
    addpath(genpath(projectRoot));
    
    % Re-verify the folders you specifically need
    addpath(genpath("module_geometry_preprocess"));
    addpath(genpath("Q3D"));
    
    fprintf('Path reset. Project Root: %s\n', projectRoot);
    
    %% ========================
    %  INITIAL POINT & BOUNDS (NOW 3 VARIABLES)
    %  ========================
    % Variables: [c_root, c_tip, alpha]
    x_0 = [1.5,  0.75,  2.0]; % Initial guess includes 2.0 degrees for alpha
    lb  = [0.8,  0.3,  -4.0]; % Lower bounds (-4 deg min alpha)
    ub  = [2.0,  0.8,  7.2]; % Upper bounds (7 deg max alpha)
        
    %% ========================
    %  OPTIMIZATION SETUP
    %  ========================
    options = optimoptions('fmincon', ...
        'Algorithm',            'sqp', ...
        'Display',              'iter', ...
        'FiniteDifferenceStepSize', 1e-2, ...
        'OutputFcn',            @my_output_fun, ...
        'MaxIterations',        100, ...
        'MaxFunctionEvaluations', 1000, ...
        'StepTolerance',        1e-4, ...
        'OptimalityTolerance',  1e-5);
        
    %% ========================
    %  RUN FMINCON
    %  ========================
    fprintf('\n--- Starting SAND optimization ---\n');
    % Notice the addition of @nonlcon_fun in the 9th argument slot
    [x_opt, fval] = fmincon(@obj_fun, x_0, [], [], [], [], lb, ub, @nonlcon_fun, options);


    %% ========================
    %  FMINCON OBJECTIVE FUNCTION
    %  ========================
    function f = obj_fun(x_bar)
        % Un-normalize variables
        c_root = x_bar(1);
        c_tip  = x_bar(2);
        alpha  = x_bar(3); % Alpha is now driven by the optimizer
        
        c           = constants();
        h_cruise    = c.altitude;
        V_fixed     = 27.78; 
        fixed_b2    = 7.5;     % 12.5m total span
        fixed_twist = 0;
    
        original_dir = pwd; 
        cleanupObj = onCleanup(@() cd(original_dir));
        
        try
            % Run the model blindly with the current fmincon alpha guess
            [LD_raw, ~, ~] = run_model([fixed_b2, c_root, c_tip, fixed_twist, V_fixed, alpha]);
            f = -LD_raw;
        catch
            f = NaN; % Backtrack if this geometry/alpha combo crashes Q3D
        end
    end

%% ========================
    %  FMINCON NONLINEAR CONSTRAINT (L = W +/- 0.5%)
    %  ========================
    function [c, ceq] = nonlcon_fun(x_bar)
        % Un-normalize variables
        c_root = x_bar(1);
        c_tip  = x_bar(2);
        alpha  = x_bar(3);
        
        c_const     = constants();
        h_cruise    = c_const.altitude;
        V_fixed     = 27.78; 
        fixed_b2    = 7.5;     
        fixed_twist = 0;
        
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
            
            % 3. The 0.5% Tolerance Band (Inequalities)
            c(1) = err - 0.005;  % Lift cannot be more than 0.5% higher
            c(2) = -err - 0.005; % Lift cannot be more than 0.5% lower
            
            % No strict equalities anymore!
            ceq = []; 
            
            % Print to console
            fprintf('    Alpha: %5.2f | Raw L_tmp: %8.2f | Target W: %8.2f | Err: %+6.2f%%\n', ...
                    alpha, L_tmp, W_target, err*100);
        catch
            % If Q3D crashes, tell fmincon the constraint is massively violated
            c = [1e6, 1e6]; 
            ceq = [];
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
            % Also print the constraint violation (how far off L=W is right now)
            fprintf('Constraint Violation : %6.4f (0 means within 0.5%% band)\n', optimValues.constrviolation);
            fprintf('======================================================\n\n');
        end
    end
end

[best_geom, best_LD] = main();