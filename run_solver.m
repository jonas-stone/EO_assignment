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


    % Calculation of Flight Conditions
    [~, rho, v, mach, mu] = calc_atmos_properties(c.altitude,v_inf,'v');
    Re = rho*v*aircraft.MAC / mu;
    aero.V = v;
    aero.rho = rho;
    aero.alt = c.altitude;
    aero.Re = Re;
    aero.M = mach;
    aero.CL = alpha;

    % Aerodynamic Solver Run
    [Res, LD, L] = Q3D_Start_mod(aircraft,aero);

    % Structural Solver Run
    structural_solver(aircraft,Res);





end