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
    %  INITIAL POINT & BOUNDS (6 VARIABLES)
    %  ========================
    % Variables: [c_root, c_tip, alpha, velocity, semi-span, twist]
    x_0 = [1.5,  0.75,  3.0, 27.78, 7.5, 0]; 
    lb  = [0.5,  0.20, -2.0, 25.00, 5.0, -5]; 
    ub  = [2.0,  0.80,  7.2, 55.55, 30.0, 0];  
    
    x_scale = [1.5, 0.75, 5.0, 30, 10, 2.0];
    
    x_0_bar = x_0 ./ x_scale;
    lb_bar  = lb  ./ x_scale;
    ub_bar  = ub  ./ x_scale;
        
    %% ========================
    %  INITIALIZE CACHE 
    %  ========================
    % This acts as our shared memory between the objective and constraint functions
    cache.x_bar = NaN(1, 6);
    cache.LD    = NaN;
    cache.L     = NaN;
    cache.W     = NaN;
    cache.stress = NaN;

    %% ========================
    %  OPTIMIZATION SETUP
    %  ========================
    options = optimoptions('fmincon', ...
        'Algorithm',            'sqp', ...
        'Display',              'iter', ...
        'FiniteDifferenceStepSize', 1e-2, ...
        'ConstraintTolerance',  1e-3, ...
        'OutputFcn',            @my_output_fun, ...
        'MaxIterations',        100, ...
        'MaxFunctionEvaluations', 1000, ...
        'StepTolerance',        1e-4, ...
        'OptimalityTolerance',  1e-5);
        
    %% ========================
    %  RUN FMINCON
    %  ========================
    fprintf('\n--- Starting CACHED COMPLEX SQP optimization ---\n');
    [x_opt, fval] = fmincon(@obj_fun, x_0_bar, [], [], [], [], lb_bar, ub_bar, @nonlcon_fun, options);
    
    % Final un-normalization for the output
    x_opt_real = x_opt .* x_scale;


    %% ========================
    %  THE CACHE ENGINE (CENTRALIZED PHYSICS)
    %  ========================
    function update_cache(x_bar)
        % 1. IF WE ALREADY SOLVED THIS EXACT GEOMETRY, DO NOTHING!
        if norm(x_bar - cache.x_bar) < 1e-10
            return; 
        end
        
        % 2. Un-normalize variables
        c_root = x_bar(1) * x_scale(1);
        c_tip  = x_bar(2) * x_scale(2);
        alpha  = x_bar(3) * x_scale(3); 
        V      = x_bar(4) * x_scale(4);
        b2     = x_bar(5) * x_scale(5);
        twist  = x_bar(6) * x_scale(6);
        
        c_const     = constants();
        h_cruise    = c_const.altitude;
        
        % Protect the directory
        original_dir = pwd; 
        cleanupObj = onCleanup(@() cd(original_dir));
        
        % Calculate Target Weight
        aircraft_tmp = calc_planform(b2, c_root, c_tip, twist);
        aero_tmp     = calc_atmos_properties(h_cruise, V, 'v', aircraft_tmp);
        [W_target, ~] = estimate_weight(aircraft_tmp, aero_tmp, V);
        
        % Update memory with the current requested state
        cache.x_bar = x_bar;
        cache.W     = W_target;
        
        % --- THE PRE-FLIGHT CAPABILITY CHECK ---
        Max_Lift_Possible = 0.5 * aero_tmp.rho * V^2 * aircraft_tmp.S * 1.5;
        
        if Max_Lift_Possible < W_target
            cache.LD = -1e6; % Objective Penalty
            cache.L  = 1e9;  % Constraint Penalty
            cache.stress = 1e9;  % Stress Penalty
            fprintf('  [CACHE] Alpha=%5.2f | Root=%5.2f | Tip=%5.2f | Sem-Span=%5.2f | Twist=%5.2f | V=%5.1f ... INSTANT FAIL (Area/Speed too low)\n', ...
                    alpha, c_root, c_tip, b2, twist, V);
            return; % Exit instantly
        end
        % ---------------------------------------
        
        fprintf('  [CACHE] Alpha=%5.2f | Root=%5.2f | Tip=%5.2f | Semi-Span=%5.2f | Twist=%5.2f | V=%5.1f ... ', ...
                alpha, c_root, c_tip, b2, twist, V);
                
        try
            % Run the heavy physics model
            [LD_raw, L_tmp, stress_tmp] = run_model([b2, c_root, c_tip, twist, V, alpha]);
            
            % Save to cache
            cache.LD = LD_raw;
            cache.L  = L_tmp;
            cache.stress = stress_tmp;
            
            err_pct = ((L_tmp - W_target) / W_target) * 100;
            fprintf('L/D = %6.2f | Err = %+6.2f%% | Stress = %6.1f MPa\n', LD_raw, err_pct, stress_tmp);
        catch
            cache.LD = -1e6; % Objective penalty
            cache.L  = 1e9;  % Constraint penalty
            cache.stress = 1e9;
            fprintf('CRASHED\n');
        end
    end

    %% ========================
    %  FMINCON OBJECTIVE FUNCTION
    %  ========================
    function f = obj_fun(x_bar)
        update_cache(x_bar); % Triggers Q3D only if needed
        f = -cache.LD;       % Grab result from memory
    end

    %% ========================
    %  FMINCON NONLINEAR CONSTRAINT (0.5% BAND)
    %  ========================
    function [c, ceq] = nonlcon_fun(x_bar)
        update_cache(x_bar); % Triggers Q3D only if needed
        
        % If Q3D crashed or failed pre-flight, output massive penalty
        if cache.L > 1e8 
            c = 1e6;
            ceq = 1e6;
        else
            % Raw fractional error from cache
            err = (cache.L - cache.W) / cache.W; 
            % The 0.5% Tolerance Band (Inequalities)
            ceq(1) = err;
                
            c_const = constants();

            % 2. STRUCTURAL CONSTRAINT
            % Carbon Fiber Compressive Limit (1200 MPa / 1.5 SF = 800 MPa)
            sigma_allow = c_const.material.Fc_y/1e6; 
            c(1) = (cache.stress / sigma_allow) - 1;
        end

    end

%% ========================
    %  CUSTOM OUTPUT FUNCTION
    %  ========================
    function stop = my_output_fun(x, optimValues, state)
        stop = false; 
        x_real = x .* x_scale;
        
       if isequal(state, 'iter')
            % Dynamically load the Ultimate Strength to calculate the margin
            c_const = constants();
            sigma_allow = c_const.material.Fc_y / 1e6; % 1200 MPa
            stress_margin = sigma_allow - cache.stress; 
            
            fprintf('\n======================================================\n');
            fprintf('ITERATION %d COMPLETED\n', optimValues.iteration);
            fprintf('Current Best L/D : %6.2f\n', -optimValues.fval); 
            fprintf('Current Geometry : Root = %5.3f | Tip = %5.3f | Semi-Span = %5.2f | Twist = %5.2f | Alpha = %5.2f | V = %5.1f\n', ...
                    x_real(1), x_real(2), x_real(5), x_real(6), x_real(3), x_real(4));
            
            fprintf('Current Stress   : %6.1f MPa (Margin: %6.1f MPa)\n', cache.stress, stress_margin);
            fprintf('Current Constr   : %6.4f (0 means valid)\n', optimValues.constrviolation);
            fprintf('======================================================\n\n');
        end
    end
end