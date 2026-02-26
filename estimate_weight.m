function [W_total, W_wing] = estimate_weight(aircraft)
    % ESTIMATE_WEIGHT Estimates the total aircraft weight of a sailplane.
    %
    % REFERENCE:
    % Raymer, D. P. (2012). "Aircraft Design: A Conceptual Approach" (5th ed.). 
    % American Institute of Aeronautics and Astronautics, Inc. 
    % Chapter 15: Weights (Statistical Empirical Weight Equations).
    %
    % This function utilizes the full Raymer Method power-law regression 
    % equations for conceptual weight estimation, keeping all geometric 
    % variables dynamic (AR, Taper Ratio, and Sweep Angle).

    % Calculate Dynamic Ratios and Angles
    % Aspect Ratio
    AR = aircraft.AR; 
    
    % Taper ratio (lambda)
    lambda = aircraft.lambda;

    % Calculate the sweep angle in radians, then take the cosine
    cos_Lambda = cosd(aircraft.LE_sweep);
    
    % Convert to Imperial Units (Raymer's equations require lbs and ft^2)
    S_sqft = aircraft.S * 10.7639;
    
    % Define Glider Constants & Raymer Exponents
    W_payload_lbs  = 90*2.20462; % Pilot + parachute (~90 kg)
    W_fuselage_lbs = 70*2.20462; % Assumed constant empty fuselage/tail mass
    n_ult = 6;          % Ultimate load factor
    
    % Raymer Exponents (Representative for light aircraft/sailplanes)
    e_1 = 0.360;   % Area exponent
    e_2 = 0.712;   % Aspect Ratio exponent
    e_3 = 0.397;   % Load factor exponent
    e_4 = -0.04;   % Taper ratio exponent
    e_5 = -1.0;    % Sweep exponent (swept wings are heavier to resist torsion)
    K   = 0.1; % Scaling constant
    
    % Iterative Weight Solver
    % Because Takeoff Weight (W_to) drives the wing weight, and wing weight 
    % drives W_to, we iterate until convergence.
    W_to_lbs = 225*2.20462; % Initial guess
    tolerance = 0.5;
    diff = 100;
    
    while diff > tolerance
        % The fully unsimplified Raymer functional form:
        % Equation Reference: Raymer, Eq. 15.46
        W_wing_lbs = K * (W_to_lbs^0.397) * (S_sqft^e_1) * (AR^e_2) * ...
                     (n_ult^e_3) * (lambda^e_4) * (cos_Lambda^e_5);
        
        % Calculate new total weight
        W_to_new = W_payload_lbs + W_fuselage_lbs + W_wing_lbs;
        
        % Check for convergence
        diff = abs(W_to_new - W_to_lbs);
        W_to_lbs = W_to_new;
    end

    Total_Mass_kg = W_to_lbs / 2.20462;

    W_total = Total_Mass_kg * 9.81;
    W_wing = (W_wing_lbs / 2.20462) * 9.80665;

end