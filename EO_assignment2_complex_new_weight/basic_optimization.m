function main()
clear all; close all; clc;
addpath("module_geometry_preprocess\");
addpath("Q3D\");

%% ========================
%  USER SETTINGS & CONSTANTS
%  ========================
run_optimization = true;
plot_interval = 1;
n_divs = 100;

% Fixed Flight and Geometric Parameters
V_fixed = 27.778;         % The strict target flight speed (m/s)
fixed_twist = 0;          % Tip twist (deg)
fixed_b2 = 7.5;          % Half-span (m) to achieve a 25m total span

% Load altitude from constants
c = constants;
h_cruise = c.altitude;

%% ========================
%  INITIAL POINT & BOUNDS (STRICTLY 2 VARIABLES)
%  ========================
% Variables: [c_root, c_tip]
x0     = [ 2,  1.5 ]; % Realistic glider chords
lb_raw = [ 0.4,  0.1 ];
ub_raw = [ 2.5,  1.5 ];

% Normalization
x0_bar = ones(1, length(x0));
lb_bar = lb_raw ./ x0;
ub_bar = ub_raw ./ x0;

% Ensure bounds aren't flipped
flipped = lb_bar > ub_bar;
temp = lb_bar(flipped);
lb_bar(flipped) = ub_bar(flipped);
ub_bar(flipped) = temp;

%% ========================
%  TRIM & EVALUATE INITIAL BASELINE
%  ========================
fprintf('--- Trimming & Evaluating initial baseline ---\n');
% We trim the baseline geometry to find its specific reference L/D and diagnostic data
[LD_init, alpha_init, trim_success, W_tot_init, W_wing_init, L_init, Cl_init] = trim_and_evaluate(x0(1), x0(2));

if ~trim_success
    error('Baseline geometry cannot support the weight at %d m/s without stalling!', V_fixed);
end

LD_ref = abs(LD_init);

fprintf('Aircraft Total Weight : %.2f N\n', W_tot_init);
fprintf('Isolated Wing Weight  : %.2f N\n', W_wing_init);
fprintf('Trimmed Alpha         : %.4f deg\n', alpha_init);
fprintf('Lift Obtained         : %.2f N\n', L_init);
fprintf('Cl Obtained           : %.4f\n', Cl_init);
fprintf('Baseline L/D          : %.2f\n', LD_init);

%% ========================
%  OPTIMIZATION SETUP
%  ========================
options = optimoptions('fmincon', ...
    'Algorithm',            'sqp', ...
    'Display',              'iter', ...
    'OutputFcn', @(x, ov, s) plot_output_fn(x, ov, s, plot_interval, n_divs, x0, fixed_b2, fixed_twist, LD_ref), ...
    'MaxIterations',        100, ...
    'MaxFunctionEvaluations', 1000, ...
    'StepTolerance',        1e-4, ...
    'OptimalityTolerance',  1e-5);

A = []; b = []; Aeq = []; beq = [];

%% ========================
%  RUN FMINCON
%  ========================
fprintf('\n--- Starting optimization ---\n');
[x_opt_bar, fval, exitflag, output] = fmincon(@obj_fun, x0_bar, ...
    A, b, Aeq, beq, lb_bar, ub_bar, [], options);

x_opt = x_opt_bar .* x0;

%% ========================
%  EVALUATE & PLOT OPTIMUM
%  ========================
fprintf('\n--- Optimum found ---\n');
% Trim the final shape one last time to get its exact parameters
[LD_opt, alpha_opt, ~, W_tot_opt, W_wing_opt, L_opt, Cl_opt] = trim_and_evaluate(x_opt(1), x_opt(2));

fprintf('Optimal c_root        = %.4f m | c_tip = %.4f m\n', x_opt(1), x_opt(2));
fprintf('Aircraft Total Weight : %.2f N\n', W_tot_opt);
fprintf('Isolated Wing Weight  : %.2f N\n', W_wing_opt);
fprintf('Required Alpha        : %.4f deg\n', alpha_opt);
fprintf('Final Lift Obtained   : %.2f N\n', L_opt);
fprintf('Final Cl Obtained     : %.4f\n', Cl_opt);
fprintf('Maximized L/D         : %.4f\n', LD_opt);
fprintf('Exit flag: %d | Iterations: %d \n', exitflag, output.iterations);

% Plot final geometry
aircraft_opt = calc_planform(fixed_b2, x_opt(1), x_opt(2), fixed_twist);
Au_opt = aircraft_opt.airfoils(1, 1:5);
Al_opt = aircraft_opt.airfoils(1, 6:10);

figure(99);
profiles_wing_3D(Au_opt, Al_opt, aircraft_opt.wing_geom, n_divs, 99);
title(sprintf('Final Optimized Wing Shape | Max L/D: %.2f', LD_opt));
drawnow;

%% ========================
%  THE CORE PHYSICS ENGINE: TRIM & EVALUATE
%  ========================
    function [best_LD, best_alpha, success, W_target, W_wing, final_L, final_Cl] = trim_and_evaluate(c_root, c_tip)
        % 1. Get the target weight
        aircraft_tmp = calc_planform(fixed_b2, c_root, c_tip, fixed_twist);
        aero_tmp = calc_atmos_properties(h_cruise, V_fixed, 'v', aircraft_tmp);
        [W_target, W_wing] = estimate_weight(aircraft_tmp, aero_tmp, V_fixed);
        
        % 2. Setup the Inner Optimizer (fminbnd)
        trim_options = optimset('Display', 'none', 'TolX', 1e-3);
        
        % Search between 0 and 8 degrees to avoid the worst extreme stalls
        [best_alpha, min_err] = fminbnd(@(a) lift_error(a, c_root, c_tip, W_target), -4, 8, trim_options);
        
        % 3. Check equilibrium
        if min_err < 0.01 
            [best_LD, final_L, ~] = run_model([fixed_b2, c_root, c_tip, fixed_twist, V_fixed, best_alpha]);
            if isnan(best_LD)
                success = false; best_LD = -inf; final_L = NaN; final_Cl = NaN;
            else
                success = true;
                
                % Calculate Lift Coefficient (Cl = L / (q * S))
                q = 0.5 * aero_tmp.rho * V_fixed^2;
                final_Cl = final_L / (q * aircraft_tmp.S);
            end
        else
            success = false; best_LD = -inf; final_L = NaN; final_Cl = NaN;
        end
    end

function err = lift_error(alpha_guess, c_root, c_tip, W_target)
        try
            % Get the raw output from the simulator
            [~, L_tmp, ~] = run_model([fixed_b2, c_root, c_tip, fixed_twist, V_fixed, alpha_guess]);
            
            % DIAGNOSTIC PRINT: See exactly what the physics engine is doing!
            fprintf('    Alpha: %5.2f | Raw L_tmp: %8.2f | Target W: %8.2f\n', alpha_guess, L_tmp, W_target);
            
            if isnan(L_tmp)
                err = 1e6 + alpha_guess^2; 
            else
                % Pure, raw comparison of Newtons to Newtons
                err = ((L_tmp - W_target) / W_target)^2; 
            end
        catch
            fprintf('    Alpha: %5.2f | CRASHED!\n', alpha_guess);
            err = 1e6 + alpha_guess^2; 
        end
    end

%% ========================
%  FMINCON OBJECTIVE FUNCTION
%  ========================
    function f = obj_fun(x_bar)
        % Un-normalize variables
        c_root = x_bar(1) * x0(1);
        c_tip  = x_bar(2) * x0(2);
        
        % Call our custom trim engine (we only need the first 3 outputs here)
        [LD_trimmed, ~, success, ~, ~, ~, ~] = trim_and_evaluate(c_root, c_tip);
        
        % Return the inverted, normalized objective to fmincon
        if success
            f = -LD_trimmed / LD_ref;
        else
            f = 10; % Penalty for geometries that stall or can't balance
        end
    end

end % <-- End of main()

%% ========================
%  PLOT OUTPUT FUNCTION
%  ========================
function stop = plot_output_fn(x_bar, optimValues, state, plot_interval, n_divs, x0_ref, fixed_b2, fixed_twist, LD_ref)
    stop = false;
    iter = optimValues.iteration;
    
    if iter == 0 || mod(iter, plot_interval) == 0
        c_root = x_bar(1) * x0_ref(1);
        c_tip  = x_bar(2) * x0_ref(2);
        
        current_LD = -optimValues.fval * LD_ref; 
        
        aircraft_iter = calc_planform(fixed_b2, c_root, c_tip, fixed_twist);
        Au = aircraft_iter.airfoils(1, 1:5);
        Al = aircraft_iter.airfoils(1, 6:10);
        
        fig_num = iter + 10; 
        figure(fig_num); 
        
        profiles_wing_3D(Au, Al, aircraft_iter.wing_geom, n_divs, fig_num);
        title(sprintf('Wing Shape — Iteration %d | Trimmed L/D: %.2f', iter, current_LD));
        drawnow;
    end
end