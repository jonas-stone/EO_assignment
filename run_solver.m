function [LD,stress_crit,L] = run_solver(x)
    c = constants();
    b2 = x(1);
    c_root = x(2);
    c_tip = x(3);
    twist_tip = x(4);
    v_inf = x(5);
    alpha = x(6);
    
    % Creating Aircraft Geometry
    aircraft = calc_planform(b2,c_root,c_tip,twist_tip);

    % Aircraft Weight Estimation
    [W_total,W_wing] = estimate_weight(aircraft);
    L = W;

    % Calculation of Flight Conditions
    [aero] = calc_atmos_properties(c.altitude,v_inf,'v',L);

    % Aerodynamic Solver Run
    [Res] = Q3D_Start_mod(aircraft,aero);

    % Structural Solver Run
    [stress_crit] = structural_solver(aircraft,Res,W_total,W_wing);





end