function [x_opt_real, fval] = main()
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
    % Variables: [c_root, c_tip, alpha, velocity, semi-span, twist]
    x_0 = [1.5,  0.75,  2.0 27.78 7.5 0]; % Initial guess includes 2.0 degrees for alpha
    lb  = [1,  0.3,  -2, 25, 5, -5]; % Lower bounds
    ub  = [2.0,  0.8,  7.2, 55.55, 30, 5]; % Upper bounds 

    x_scale = [1.5, 0.75, 5.0, 30, 10, 2.0];
    x_0_bar =  x_0./ x_scale;
    lb_bar = lb ./ x_scale;
    ub_bar = ub ./ x_scale;
        
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
    fprintf('\n--- Starting COMPLEX SQP optimization ---\n');
    % Notice the addition of @nonlcon_fun in the 9th argument slot
    [x_opt, fval] = fmincon(@obj_fun, x_0_bar, [], [], [], [], lb_bar, ub_bar, @nonlcon_fun, options);

    % Final un-normalization for the output
    x_opt_real = x_opt .* x_scale;

%% ========================
    %  FMINCON OBJECTIVE FUNCTION
    %  ========================
    function f = obj_fun(x_bar)
        % Un-normalize variables using multiplication (*)
        c_root = x_bar(1) * x_scale(1);
        c_tip  = x_bar(2) * x_scale(2);
        alpha  = x_bar(3) * x_scale(3); 
        V      = x_bar(4) * x_scale(4);
        b2     = x_bar(5) * x_scale(5);
        twist  = x_bar(6) * x_scale(6);
        
        fprintf('  [OBJ] Testing Alpha = %5.2f | Root = %5.2f | Tip = %5.2f | Span = %5.2f | Twist = %5.2f | V = %5.1f ... ', ...
                alpha, c_root, c_tip, b2, twist, V);
                
        c           = constants();
        h_cruise    = c.altitude;
        original_dir = pwd; 
        cleanupObj = onCleanup(@() cd(original_dir));
        
        % --- THE PRE-FLIGHT CAPABILITY CHECK ---
        % Calculate Wing Area
        S = b2 * (c_root + c_tip); 
        
        % Calculate approximate air density at 1500m
        rho = 1.0581; 
        
        % Calculate the absolute Maximum Lift this geometry could ever produce
        % assuming a highly optimistic Max CL of 1.5
        Max_Lift_Possible = 0.5 * rho * V^2 * S * 1.5;
        
        % We need the rough target weight to compare it to
        aircraft_tmp = calc_planform(b2, c_root, c_tip, twist);
        aero_tmp     = calc_atmos_properties(h_cruise, V, 'v', aircraft_tmp);
        [W_target, ~] = estimate_weight(aircraft_tmp, aero_tmp, V);
        
        % If the maximum possible lift is LESS than the weight, this combo is impossible!
        if Max_Lift_Possible < W_target
            f = 1e6;
            fprintf('INSTANT FAIL (Area/Speed too low)\n');
            return; % Exit the function instantly without running Q3D!
        end
        % ---------------------------------------
        
        try
            % If it passes the check, run the heavy physics model
            [LD_raw, ~, ~] = run_model([b2, c_root, c_tip, twist, V, alpha]);
            f = -LD_raw;
            fprintf('L/D = %6.2f\n', LD_raw);
        catch
            f = 1e6; 
            fprintf('CRASHED\n');
        end
    end

    %% ========================
    %  FMINCON NONLINEAR CONSTRAINT (0.5% BAND)
    %  ========================
    function [c, ceq] = nonlcon_fun(x_bar)
        % Un-normalize variables using multiplication (*)
        c_root = x_bar(1) * x_scale(1);
        c_tip  = x_bar(2) * x_scale(2);
        alpha  = x_bar(3) * x_scale(3); 
        V      = x_bar(4) * x_scale(4);
        b2     = x_bar(5) * x_scale(5);
        twist  = x_bar(6) * x_scale(6);
        
        c_const     = constants();
        h_cruise    = c_const.altitude;
        
        original_dir = pwd; 
        cleanupObj = onCleanup(@() cd(original_dir));
        
        % 1. Calculate Target Weight for this geometry
        aircraft = calc_planform(b2, c_root, c_tip, twist);
        aero     = calc_atmos_properties(h_cruise, V, 'v', aircraft);
        [W_target, ~] = estimate_weight(aircraft, aero, V);
        
        % --- THE PRE-FLIGHT CAPABILITY CHECK ---
        % Calculate Wing Area
        S = b2 * (c_root + c_tip); 
        
        % Calculate approximate air density at 1500m
        rho = 1.0581; 
        
        % Calculate the absolute Maximum Lift this geometry could ever produce
        % assuming a highly optimistic Max CL of 1.5
        Max_Lift_Possible = 0.5 * rho * V^2 * S * 1.5;
        
        % If the maximum possible lift is LESS than the weight, this combo is impossible!
        if Max_Lift_Possible < W_target
            c = [1e6, 1e6]; % Massive constraint violation
            ceq = [];
            fprintf('    Alpha: %5.2f | Raw L_tmp:  INSTANT FAIL (Area/Speed too low) | Target W: %8.2f\n', alpha, W_target);
            return; % Exit the function instantly without running Q3D!
        end
        % ---------------------------------------
        
        % 2. Calculate Actual Lift at the guessed alpha
        try
            [~, L_tmp, ~] = run_model([b2, c_root, c_tip, twist, V, alpha]);
            
            % Raw fractional error
            err = (L_tmp - W_target) / W_target; 
            
            % The 0.5% Tolerance Band (Inequalities)
            c(1) = err - 0.005;  % Lift cannot be more than 0.5% higher
            c(2) = -err - 0.005; % Lift cannot be more than 0.5% lower
            
            ceq = []; % Leave equality blank
            
            fprintf('    Alpha: %5.2f | Raw L_tmp: %8.2f | Target W: %8.2f | Err: %+6.2f%%\n', ...
                    alpha, L_tmp, W_target, err*100);
        catch
            c = [1e6, 1e6]; % Massive violation penalty assigned to 'c'
            ceq = [];
            fprintf('    Alpha: %5.2f | Raw L_tmp:  CRASHED | Target W: %8.2f\n', alpha, W_target);
        end
    end

    %% ========================
    %  CUSTOM OUTPUT FUNCTION
    %  ========================
    function stop = my_output_fun(x, optimValues, state)
        stop = false; 
        % Un-normalize for readable printing
        x_real = x .* x_scale;
        
       if isequal(state, 'iter')
            fprintf('\n======================================================\n');
            fprintf('ITERATION %d COMPLETED\n', optimValues.iteration);
            fprintf('Current Best L/D : %6.2f\n', -optimValues.fval); 
            fprintf('Current Geometry : Root = %5.3f | Tip = %5.3f | Span = %5.2f | Twist = %5.2f | Alpha = %5.2f | V = %5.1f\n', ...
                    x_real(1), x_real(2), x_real(5), x_real(6), x_real(3), x_real(4));
            fprintf('Current Constr   : %6.4f (0 means within 0.5%% band)\n', optimValues.constrviolation);
            fprintf('======================================================\n\n');
        end
    end
end

[best_geom, best_LD] = main();