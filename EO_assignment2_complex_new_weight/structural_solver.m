function [stress_crit] = structural_solver(aircraft, Res, W_total, W_wing_total, W_water_total)
    % 1. Define Ultimate Load Condition
    n_ult = 7.95; % JAR-22 Utility Category Ultimate Load Factor
    
    % Apply load factor to all forces
    L_half_req = (W_total * n_ult) / 2; 
    W_wing_half = (W_wing_total * n_ult) / 2;
    W_water_half = (W_water_total * n_ult) / 2; 

    % 2. Extract Q3D Spanwise Data
    Y = Res.Section.Y(:)';   
    CL = Res.Section.Cl(:)'; 

    % 3. Calculate Local Chord along the span
    y_root = aircraft.wing_geom(1,2);
    y_tip  = aircraft.wing_geom(2,2);
    c_root = aircraft.wing_geom(1,4);
    c_tip  = aircraft.wing_geom(2,4);
    
    c_y = interp1([y_root, y_tip], [c_root, c_tip], Y);

    % 4. Aerodynamic Lift Distribution (Upward Force)
    L_shape = CL .* c_y;
    Area_shape = trapz(Y, L_shape);
    l_y = L_shape * (L_half_req / Area_shape); % Scaled Lift (N/m)

    % 5. Downward Inertial Relief (Wing Mass + Water Ballast)
    Area_chord = trapz(Y, c_y);
    w_wing_y = c_y * (W_wing_half / Area_chord); 
    w_water_y = c_y * (W_water_half / Area_chord); 
    w_total_down_y = w_wing_y + w_water_y; 

    % 6. Net Force and Root Bending Moment Integration
    net_force_y = l_y - w_total_down_y;
    
    % Ensure Y is in ascending order (Root to Tip) before flipping
    % Flip arrays to integrate from the free TIP inward to the fixed ROOT
    Y_flipped = flip(Y);
    force_flipped = flip(net_force_y);
    
    % Step A: Integrate Load to get Shear Force (V)
    % Note: Since we are integrating backwards in space (dx is negative), 
    % we take the absolute value of the diff(Y) inside cumtrapz to keep signs correct,
    % or simply use the flipped Y array which MATLAB handles as negative steps.
    % To keep it mathematically simple and positive:
    dy_flipped = abs(diff(Y_flipped));
    V_shear = zeros(1, length(Y_flipped));
    for i = 2:length(Y_flipped)
        % Trapezoidal integration step
        V_shear(i) = V_shear(i-1) + 0.5 * (force_flipped(i) + force_flipped(i-1)) * dy_flipped(i-1);
    end
    
    % Step B: Integrate Shear Force to get Bending Moment (M)
    M_dist = zeros(1, length(Y_flipped));
    for i = 2:length(Y_flipped)
        % Trapezoidal integration step
        M_dist(i) = M_dist(i-1) + 0.5 * (V_shear(i) + V_shear(i-1)) * dy_flipped(i-1);
    end
    
    % The Root Bending Moment is the final accumulated value
    M_root = M_dist(end);

    % 7. Define the Wing Root Structural Cross-Section
    tc_ratio = 0.1438; 
    t_max = c_root * tc_ratio; 
    
    % REALISTIC SPAR DIMENSIONS
    box_width = 0.2 * c_root; 
    box_height = t_max*0.9; 
    
    % HARDCODED SPAR DIMENSIONS
    t_web = 0.003; % 3mm shear webs
    t_cap_root = 0.008; % carbon spar cap at the root
    
    % 8. Calculate Second Moment of Area (I_xx)
    B_outer = box_width;
    H_outer = box_height;
    B_inner = box_width - 2*t_web;
    H_inner = box_height - 2*t_cap_root;
    
    I_xx = (1/12) * B_outer * H_outer^3 - (1/12) * B_inner * H_inner^3;
    
    % 9. Calculate Critical Bending Stress
    y_max = H_outer / 2; 
    stress_crit_Pa = (M_root * y_max) / I_xx; 
    stress_crit = stress_crit_Pa / 1e6; % Convert to MPa
end