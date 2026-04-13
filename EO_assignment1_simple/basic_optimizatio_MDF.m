function main()
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
    %  1. SETTINGS & SWEEP RANGE
    %  ========================
    % Fixed Parameters (Matches your 12.5m span glider)


    %% ========================
    %  INITIAL POINT & BOUNDS (STRICTLY 2 VARIABLES)
    %  ========================
    % Variables: [c_root, c_tip]
    x_0 = [1.5, 0.75]; % Chord within the possible geometries
    lb = [0.8, 0.2];
    ub = [2, 0.8];
        
%% ========================
%  OPTIMIZATION SETUP
%  ========================
options = optimoptions('fmincon', ...
    'Algorithm',            'sqp', ...
    'Display',              'iter', ...
    'FiniteDifferenceStepSize', 1e-2, ...
    'OutputFcn',           @my_output_fun, ...
    'MaxIterations',        100, ...
    'MaxFunctionEvaluations', 1000, ...
    'StepTolerance',        1e-4, ...
    'OptimalityTolerance',  1e-5);

%% ========================
%  RUN FMINCON
%  ========================
fprintf('\n--- Starting optimization ---\n');
[x_opt, fval] = fmincon(@obj_fun, x_0, [], [], [], [], lb, ub, [], options);

%% ========================
%  INTERNAL PHYSICS ENGINE
%  ========================
function [best_LD, best_alpha, success] = trim_engine(c_root, c_tip, b2, twist, V, h)
    % 1. Get target weight for this shape
    aircraft = calc_planform(b2, c_root, c_tip, twist);
    aero     = calc_atmos_properties(h, V, 'v', aircraft);
    [W_target, ~] = estimate_weight(aircraft, aero, V);
    
    % 2. Solve for Alpha (Trim)
    % We allow -4 to 8 deg for the sweep to ensure we find a solution
    opts = optimset('Display', 'none', 'TolX', 1e-3);
    [best_alpha, min_err] = fminbnd(@(a) lift_err(a, b2, c_root, c_tip, twist, V, W_target), -6, 10, opts);
    
    % 3. Extract final L/D
    if min_err < 0.0001 
        [best_LD, ~, ~] = run_model([b2, c_root, c_tip, twist, V, best_alpha]);
        if isnan(best_LD), success = false; else, success = true; end
    else
        success = false;
        best_LD = 0;
    end
end

function err = lift_err(a, b2, cr, ct, tw, V, W_t)
    % 1. Remember exactly where we are before we run the physics engine
    original_dir = pwd; 
    cleanupObj = onCleanup(@() cd(original_dir));
    try
        % Run the model first (silently)
        [~, L_tmp, ~] = run_model([b2, cr, ct, tw, V, a]);
        
        % 2. If it survives, print your exact desired line!
        fprintf('    Alpha: %5.2f | Raw L_tmp: %8.2f | Target W: %8.2f\n', a, L_tmp, W_t);
        
        % Calculate error
        if isnan(L_tmp)
            err = 1e6; 
        else
            err = ((L_tmp - W_t)/W_t)^2; 
        end
        
    catch ME
        % 3. If it crashes, print a matching line to keep the log uniform
        fprintf('    Alpha: %5.2f | Raw L_tmp:  CRASHED | Target W: %8.2f\n', a, W_t);
       
        % Assign a massive error penalty so fminbnd knows this alpha is invalid
        err = 1e6;
    end
end

%% ========================
%  FMINCON OBJECTIVE FUNCTION
%  ========================
    function f = obj_fun(x_bar)
        % Un-normalize variables
        c_root = x_bar(1);
        c_tip  = x_bar(2);
        c           = constants();
        h_cruise    = c.altitude;
        V_fixed     = 27.78; 
        fixed_b2    = 7.5;     % 12.5m total span
        fixed_twist = 0;
    
        % Call our custom trim engine
        [LD_trimmed, ~, success] = trim_engine(c_root, c_tip, fixed_b2, fixed_twist, V_fixed, h_cruise);
        
        % Return the inverted, normalized objective to fmincon
        if success
            f = -LD_trimmed;
        else
            f = NaN; 
        end
    end

%% ========================
%  CUSTOM OUTPUT FUNCTION
%  ========================
function stop = my_output_fun(x, optimValues, state)
    stop = false; % This must be here, it tells fmincon not to abort
    
    % Only print when fmincon officially completes an iteration step
    if isequal(state, 'iter')
        fprintf('\n======================================================\n');
        fprintf('ITERATION %d COMPLETED\n', optimValues.iteration);
        fprintf('Current Best L/D : %6.2f\n', -optimValues.fval); % Un-invert the objective
        fprintf('Current Geometry : c_root = %5.3f | c_tip = %5.3f\n', x(1), x(2));
        fprintf('======================================================\n\n');
    end
end

end

