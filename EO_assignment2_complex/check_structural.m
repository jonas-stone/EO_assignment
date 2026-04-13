% --- STANDALONE STRUCTURAL TEST SCRIPT ---
clear; clc; close all;

% 1. Define a basic wing (Span=7.5, Root=1.5, Tip=0.75, Twist=0, V=27.8, Alpha=2)
b2 = 7.5; c_root = 1.5; c_tip = 0.75; twist = 0; V = 27.8; alpha = 4.0;
h_cruise = 1500;

% 2. Build the basic structs
aircraft = calc_planform(b2, c_root, c_tip, twist);
aero = calc_atmos_properties(h_cruise, V, 'v', aircraft);
aero.Alpha = alpha; % Ensure alpha is passed in

% 3. Run ONLY the Q3D solver once
fprintf('Running Q3D Physics...\n');
[Res, LD, L] = Q3D_Start_mod(aircraft, aero);

% 4. Estimate Weights (We need these for the structure!)
fprintf('Estimating Weights...\n');
[W_total, W_wing_total, W_water_total] = estimate_weight(aircraft, aero, V);


% 5. Run the Structural Solver
fprintf('Calculating Root Bending Stress...\n');
stress_crit = structural_solver(aircraft, Res, W_total, W_wing_total, W_water_total);

% 6. Print the Final Results
fprintf('\n======================================\n');
fprintf('          TEST RESULTS\n');
fprintf('======================================\n');
fprintf('Aircraft Weight : %8.2f N\n', W_total);
fprintf('1G Lift Made    : %8.2f N\n', L);
fprintf('Critical Stress : %8.2f MPa\n', stress_crit);
fprintf('Allowed Stress  : %8.2f MPa (assuming Carbon Yield 276/1.5)\n', 184.00);

if stress_crit > 184
    fprintf('Result          : WING SNAPPED!\n');
else
    fprintf('Result          : WING SURVIVED!\n');
end
fprintf('======================================\n');