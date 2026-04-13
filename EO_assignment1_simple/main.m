function main()   % <-- wrap everything in a function
clear all; close all; clc;

addpath("module_geometry_preprocess\");
addpath("Q3D\");

%% ========================
%  USER SETTINGS
%  ========================
run_optimization = true;
plot_interval = 5;
n_divs = 100;

%% ========================
%  INITIAL POINT & BOUNDS
%  ========================
% b2, c_root, c_tip, twist_tip, v_inf, alpha
x0     = [ 7.5,  0.9,   0.35,   -2,       25,    3  ];
lb_raw = [ 1,  0.5,   0.1,   -10,       5,    0  ];
ub_raw = [ 30, 5,   5,    10, 100,    8  ];

x0_bar = ones(1, length(x0));
lb_bar = lb_raw ./ x0;
ub_bar = ub_raw ./ x0;

flipped = lb_bar > ub_bar;
temp = lb_bar(flipped);
lb_bar(flipped) = ub_bar(flipped);
ub_bar(flipped) = temp;

%% ========================
%  EVALUATE & PLOT INITIAL
%  ========================
fprintf('--- Evaluating initial point ---\n');
[LD_init, L_init, stress_init] = run_model(x0);
fprintf('Initial L/D = %.4f | L = %.2f N | stress_crit = %.2f MPa\n', ...
    LD_init, L_init, stress_init);

% 1. Create the aircraft geometry
aircraft_init = calc_planform(x0(1), x0(2), x0(3), x0(4));

c = constants;
v_inf = x0(5);
h_cruise = c.altitude;

% 2. Calculate Aero Properties (Requires the aircraft struct for MAC)
aero_init = calc_atmos_properties(h_cruise, v_inf, 'v', aircraft_init);

% 3. Estimate Weight (Pass the calculated aero and V_inf)
[W_total_init, ~] = estimate_weight(aircraft_init, aero_init, aero_init.V);

% Normalization references
LD_ref     = abs(LD_init);
L_ref      = W_total_init;
sigma_yield = 276e6;
SF = 1.5;
sigma_allow = sigma_yield / SF;
stress_ref  = sigma_allow;

Au = aircraft_init.airfoils(1, 1:5);
Al = aircraft_init.airfoils(1, 6:10);
profiles_wing_3D(Au, Al, aircraft_init.wing_geom, n_divs, 1);
title('Initial Wing Shape');
drawnow;


%% ========================
%  INITIALIZE CACHE (before fmincon!)
%  ========================
cache.x      = NaN(size(x0));
cache.LD     = NaN;
cache.L      = NaN;
cache.stress = NaN;
cache.W      = NaN;

%% ========================
%  OPTIMIZATION SETUP
%  ========================
options = optimoptions('fmincon', ...
    'Algorithm',            'sqp', ...
    'Display',              'iter', ...
    'OutputFcn',            @(x, ov, s) plot_output_fn(x, ov, s, plot_interval, n_divs, x0), ...
    'MaxIterations',        200, ...
    'MaxFunctionEvaluations', 2000, ...
    'StepTolerance',        1e-6, ...
    'OptimalityTolerance',  1e-6);

A = []; b = []; Aeq = []; beq = [];

%% ========================
%  RUN FMINCON
%  ========================
fprintf('\n--- Starting optimization ---\n');
[x_opt_bar, fval, exitflag, output] = fmincon(@obj_fun, x0_bar, ...
    A, b, Aeq, beq, lb_bar, ub_bar, @nonlcon, options);

x_opt = x_opt_bar .* x0;

%% ========================
%  EVALUATE & PLOT OPTIMUM
%  ========================
[LD_opt, L_opt, stress_opt] = run_model(x_opt);
fprintf('\n--- Optimum found ---\n');
fprintf('Optimal L/D = %.4f | L = %.2f N | stress_crit = %.2f MPa\n', ...
    LD_opt, L_opt, stress_opt);
fprintf('x_opt = [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n', x_opt);
fprintf('Exit flag: %d | Iterations: %d | f-evals: %d\n', ...
    exitflag, output.iterations, output.funcCount);

aircraft_opt = calc_planform(x_opt(1), x_opt(2), x_opt(3), x_opt(4));
Au_opt = aircraft_opt.airfoils(1, 1:5);
Al_opt = aircraft_opt.airfoils(1, 6:10);
profiles_wing_3D(Au_opt, Al_opt, aircraft_opt.wing_geom, n_divs, 1);
title('Optimized Wing Shape');
drawnow;

%% ========================
%  NESTED FUNCTIONS (share parent workspace)
%  ========================
    function update_cache(x_bar)
        x_real = x_bar .* x0;
        if isequal(x_bar, cache.x)
            return;
        end
        
        [tmpLD, tmpL, tmpStress] = run_model(x_real);
        cache.LD     = tmpLD;
        cache.L      = tmpL;
        cache.stress = tmpStress;
        
        % 1. Get new geometry
        aircraft_tmp = calc_planform(x_real(1), x_real(2), x_real(3), x_real(4));
        
        % 2. Extract the current velocity from the active design variables
        v_inf_tmp = x_real(5); 
        
        % 3. Recalculate aero properties for the new geometry & speed
        % (h_cruise is accessible because it's defined in the main workspace)
        aero_tmp = calc_atmos_properties(h_cruise, v_inf_tmp, 'v', aircraft_tmp);
        
        % 4. Estimate weight using ALL THREE required arguments
        [tmpW, ~] = estimate_weight(aircraft_tmp, aero_tmp, v_inf_tmp);
        
        cache.W = tmpW;
        cache.x = x_bar;
    end

    function f = obj_fun(x_bar)
        update_cache(x_bar);
        f = -cache.LD / LD_ref;
    end

    function [c, ceq] = nonlcon(x_bar)
        update_cache(x_bar);
        c(1) = (cache.stress - sigma_allow) / stress_ref;
        ceq(1) = (cache.L - cache.W) / L_ref;
    end

end  % end of function main()

%% ========================
%  STANDALONE LOCAL FUNCTION (doesn't need workspace access)
%  ========================
function stop = plot_output_fn(x_bar, optimValues, state, plot_interval, n_divs, x0_ref)
    stop = false;
    iter = optimValues.iteration;
    if iter == 0 || mod(iter, plot_interval) == 0
        x_real = x_bar .* x0_ref;
        aircraft_iter = calc_planform(x_real(1), x_real(2), x_real(3), x_real(4));
        Au = aircraft_iter.airfoils(1, 1:5);
        Al = aircraft_iter.airfoils(1, 6:10);
        profiles_wing_3D(Au, Al, aircraft_iter.wing_geom, n_divs, 1);
        title(sprintf('Wing Shape — Iteration %d', iter));
        drawnow;
    end
end