function [LD, L, stress_crit] = run_model(x)
    c = constants();
    
    % Your original 6 design variables
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
    
    % Aircraft Weight Estimation (W_total_empty = Wing + Fuselage + Pilot)
    [W_total_empty, W_wing] = estimate_weight(aircraft, aero, v_inf);
    
    % --- WATER BALLAST LOGIC ---
    % MTOW based on EASA TCDS Nimbus-4M limit (800 kg)
    W_MTOW = 800;
    
    % Available water ballast is whatever weight is left over
%   W_water_total = W_MTOW - W_total_empty; 
%    if W_water_total < 0
%       W_water_total = 0; % Prevent negative water if the wing gets too heavy
%       W_MTOW = W_total_empty; 
%   end
    
    % Aerodynamic Solver Run
    [Res, LD, L] = Q3D_Start_mod(aircraft, aero);
    
   % if ~isnan(LD)
        % Structural Solver Run 
        % We pass W_MTOW because the worst-case structural load is at max weight
  %      [stress_crit] = structural_solver(aircraft, Res, W_MTOW, W_wing, W_water_total);
  %  else
     stress_crit = NaN;
   % end
end