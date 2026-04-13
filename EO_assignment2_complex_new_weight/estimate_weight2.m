function [W_total, W_wing, W_water] = estimate_weight2(aircraft, aero, V_inf)
    % ESTIMATE_WEIGHT Estimates the total aircraft weight of a sailplane.
    %
    % REFERENCE:
    % Raymer, D. P. (2012). "Aircraft Design: A Conceptual Approach" (5th ed.). 
    % American Institute of Aeronautics and Astronautics, Inc. 
    % Chapter 15: Weights (Statistical Empirical Weight Equations).

    % Aspect Ratio
    AR = aircraft.AR; 
    
    % Taper ratio (lambda)
    lambda = aircraft.lambda;

    % Calculate the sweep angle in radians, then take the cosine
    cos_Lambda = cosd(aircraft.LE_sweep);
    
    % Convert to Imperial Units (Raymer's equations require lbs and ft^2)
    S_sqft = aircraft.S * 10.7639;

    % Density [kg/m^3]
    rho = aero.rho;

    % Dynamic pressure at cruise lb/ft^2
    q = (0.5 * rho * V_inf^2) * 0.0208854;
    
    % Airfoil thickness to chord ratio
    t_c = 0.14;
    
    % Ultimate load factor 
    n_ult = 7.95; 

    % Maximum take-off weight (in lbs)
    W_max = 550*2.20462; 

    % Payload and weight
    M_payload  = 90; % Pilot + parachute (~90 kg)
    M_non_lifting = 240; % Assumed constant empty fuselage/tail mass
    
    % Raymer Exponents
    k = 0.036;
    e_1 = 0.758;
    e_2 = 0.6;
    e_3 = 0.006;
    e_4 = 0.04;
    e_5 = -0.3;
    e_6 = 0.49;

    % Mass of the wing [lb]
    W_wing_imperial = k*((S_sqft^e_1)*((AR/(cos_Lambda^2))^e_2)*(q^e_3)*(lambda^e_4)*(((100*t_c)/cos_Lambda)^e_5)*(n_ult*W_max)^e_6);
    W_wing = (W_wing_imperial/2.20462)*9.81;
    
    % --- NEW: Water Ballast Assumption ---
    % Set to 0 N for "dry" flight, or add weight here (e.g., 50 * 9.81 for 50 liters)
    W_water = 0*9.81;

    % Total aircraft weight [N]
    W_total = W_wing + (M_payload + M_non_lifting)*9.81 + W_water;

end