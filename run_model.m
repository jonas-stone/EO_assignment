function [LD,L,stress_crit] = run_model(x)
    c = constants();
    b2 = x(1);
    c_root = x(2);
    c_tip = x(3);
    twist_tip = x(4);
    v_inf = x(5);
    alpha = x(6);
    
    % Creating Aircraft Geometry
    aircraft = calc_planform(b2,c_root,c_tip,twist_tip);

    % Calculation of Flight Conditions
    [aero] = calc_atmos_properties(c.altitude,v_inf,'v',aircraft);
    aero.Alpha = alpha;

    % Aircraft Weight Estimation
    [W_total,W_wing] = estimate_weight(aircraft, aero, v_inf);

    % Aerodynamic Solver Run
    [Res, LD, L] = Q3D_Start_mod(aircraft,aero);
    
    if ~isnan(LD)
        % Structural Solver Run
        [stress_crit] = structural_solver(aircraft,Res,W_total,W_wing);
    else
        stress_crit = NaN;
    end

end