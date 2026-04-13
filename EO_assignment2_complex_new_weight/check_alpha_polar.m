function check_alpha_polar()
    close all; clc;
    
    % --- Setup Paths ---
    projectRoot = pwd; 
    addpath(genpath(projectRoot));
    addpath(genpath("module_geometry_preprocess"));
    addpath(genpath("Q3D"));
    
    % --- The Problematic Geometry ---
    c_root = 0.8;
    c_tip  = 0.41;
    b2     = 7.5;      
    twist  = 0;
    V      = 35;    
    
    % --- EXACT DATA EXTRACTION ---
    % Ask your own engine for the parameters so nothing is assumed
    c_const  = constants();
    aircraft = calc_planform(b2, c_root, c_tip, twist);
    aero     = calc_atmos_properties(c_const.altitude, V, 'v', aircraft);
    
    % Extract density and area from your structs 
    % (Note: If your struct fields are named slightly differently, e.g., aero.density, change them here!)
   rho = aero.rho;

   S = aircraft.S;
    
    q = 0.5 * rho * V^2; % Dynamic pressure
    
    % --- Sweep Parameters ---
    alpha_vec = 5:0.2:8;
    CL_vec = zeros(size(alpha_vec));
    L_vec  = zeros(size(alpha_vec));
    converged_flags = true(size(alpha_vec));
    
    fprintf('Starting Alpha Sweep for Root=%.4f, Tip=%.4f\n', c_root, c_tip);
    fprintf('Area = %.2f m^2 | Density = %.4f kg/m^3\n', S, rho);
    fprintf('------------------------------------------------------------\n');
    fprintf(' Alpha (deg) |   Lift (N)  |    CL    | Status \n');
    fprintf('------------------------------------------------------------\n');
    
    for i = 1:length(alpha_vec)
        a = alpha_vec(i);
        
        % Clear the last warning so we can catch a fresh one
        lastwarn('');
        
        tic; % Start timer
        try
            % Call YOUR exact run_model (included at the bottom)
            [~, L_tmp, ~] = run_model([b2, c_root, c_tip, twist, V, a]);
            
            % Check if Q3D generated a warning
            [warnMsg, ~] = lastwarn;
            if contains(warnMsg, 'did not converge') || isnan(L_tmp)
                converged_flags(i) = false;
            end
            
            % Calculate CL
            L_vec(i) = L_tmp;
            CL_vec(i) = L_tmp / (q * S);
            elapsed = toc;
            
            % Print formatted results
            if converged_flags(i)
                fprintf(' %11.2f | %11.2f | %8.4f | OK (%.1fs)\n', a, L_vec(i), CL_vec(i), elapsed);
            else
                fprintf(' %11.2f | %11.2f | %8.4f | FAILED TO CONVERGE (%.1fs)\n', a, L_vec(i), CL_vec(i), elapsed);
            end
            
        catch
            converged_flags(i) = false;
            L_vec(i) = NaN;
            CL_vec(i) = NaN;
            elapsed = toc;
            fprintf(' %11.2f | %11s | %8s | CRASHED (%.1fs)\n', a, 'NaN', 'NaN', elapsed);
        end
    end
    
    % --- Visualization ---
    figure('Name', 'Lift Coefficient Polar', 'Color', 'w');
    hold on; grid on;
    
    % Plot valid points (Blue Line)
    valid_idx = converged_flags & ~isnan(CL_vec);
    plot(alpha_vec(valid_idx), CL_vec(valid_idx), 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'DisplayName', 'Healthy Flow');
    
    % Plot unconverged points (Red Triangles)
    invalid_idx = ~converged_flags & ~isnan(CL_vec);
    if any(invalid_idx)
        plot(alpha_vec(invalid_idx), CL_vec(invalid_idx), 'r^', 'MarkerSize', 9, 'MarkerFaceColor', 'r', 'LineWidth', 1.5, 'DisplayName', 'Q3D Convergence Warning');
    end
    
    title(sprintf('C_L vs \\alpha Polar\n(Root: %.4f m, Tip: %.4f m)', c_root, c_tip));
    xlabel('Angle of Attack, \alpha (deg)');
    ylabel('Lift Coefficient, C_L');
    legend('Location', 'northwest');
    xlim([-5 11]);
    
    %% ========================
    %  YOUR EXACT RUN_MODEL
    %  ========================
    function [LD, L, stress_crit] = run_model(x)
        c = constants();
        
        b2 = x(1);
        c_root = x(2);
        c_tip = x(3);
        twist_tip = x(4);
        v_inf = x(5);
        alpha = x(6);
        
        % Creating Aircraft Geometry
        aircraft = calc_planform(b2, c_root, c_tip, twist_tip);
        
        % Calculation of Flight Conditions
        [aero] = calc_atmos_properties(c.altitude, v_inf, 'v', aircraft);
        aero.Alpha = alpha;
        
        % Aircraft Weight Estimation
        [W_total_empty, W_wing] = estimate_weight(aircraft, aero, v_inf);
        
        % Aerodynamic Solver Run
        [Res, LD, L] = Q3D_Start_mod(aircraft, aero);
        
        stress_crit = NaN;
    end
end