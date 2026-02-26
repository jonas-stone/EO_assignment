function [stress_crit] = structural_solver(aircraft, Res, W_total, W_wing_total)
    % STRUCTURAL_SOLVER Calculates critical root bending stress by numerically 
    % integrating the Q3D spanwise lift distribution and applying inertial relief.

    % We are analyzing one half of the wing
    L_half_req = W_total / 2; 
    W_wing_half = W_wing_total / 2;

    % 2. Extract Q3D Spanwise Data
    Y = Res.Section.Y;   % Spanwise coordinates (m)
    CL = Res.Section.Cl; % Local lift coefficient distribution

    % Ensure Y and CL are row vectors for consistent matrix math
    Y = Y(:)';
    CL = CL(:)';

    % 3. Calculate Local Chord along the span (c_y)
    y_root = aircraft.wing_geom(1,2);
    y_tip  = aircraft.wing_geom(2,2);
    c_root = aircraft.wing_geom(1,4);
    c_tip  = aircraft.wing_geom(2,4);

    % Linearly interpolate the chord at the exact Y stations Q3D used
    c_y = interp1([y_root, y_tip], [c_root, c_tip], Y);

    % 4. Aerodynamic Lift Distribution (Upward Force)
    % Lift shape is proportional to local CL * local chord.
    L_shape = CL .* c_y;
    
    % Integrate the shape to find the area under the curve
    Area_shape = trapz(Y, L_shape);
    
    % Scale the distribution so the total integrated lift perfectly equals L_half_req
    l_y = L_shape * (L_half_req / Area_shape); % Lift per meter (N/m)

    % 5. Wing Weight Distribution / Inertial Relief (Downward Force)
    % The structural mass is roughly proportional to the local chord area.
    Area_chord = trapz(Y, c_y);
    w_y = c_y * (W_wing_half / Area_chord); % Weight per meter (N/m)

    % 6. Net Force and Root Bending Moment Integration
    % Net vertical force per meter (Upward Lift - Downward Weight)
    net_force_y = l_y - w_y;

    % Numerically integrate (Force * moment_arm) to get Root Bending Moment.
    % For the root, the moment arm to any station Y is simply Y.
    integrand = net_force_y .* Y;
    M_root = trapz(Y, integrand); % Final Root Bending Moment in N*m

    % 7. Define the Wing Root Structural Cross-Section (Simplified Box Beam)
    tc_ratio = 0.13; % Assumed 13% airfoil thickness
    t_max = c_root * tc_ratio; 
    
    box_width = 0.40 * c_root; % Box occupies 40% of the chord
    box_height = t_max; 
    
    t_cap = 0.005; % Spar cap thickness (5 mm)
    t_web = 0.003; % Spar web thickness (3 mm)
    
    % 8. Calculate Second Moment of Area (I_xx)
    B_outer = box_width;
    H_outer = box_height;
    B_inner = box_width - 2*t_web;
    H_inner = box_height - 2*t_cap;
    
    I_xx = (1/12) * B_outer * H_outer^3 - (1/12) * B_inner * H_inner^3;
    
    % 9. Calculate Critical Bending Stress (Navier's Equation)
    y_max = H_outer / 2; 
    stress_crit_Pa = (M_root * y_max) / I_xx; 
    
    % Convert to MegaPascals (MPa)
    stress_crit = stress_crit_Pa / 1e6; 
    
end