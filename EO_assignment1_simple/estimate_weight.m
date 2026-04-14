function [W_total, W_wing] = estimate_weight(aircraft, aero, V_inf)
    % HYBRID COMPONENT BUILD-UP WEIGHT ESTIMATION
    
    % 1. Aircraft Geometry Extraction
    b2 = aircraft.b2; % Semi-span
    c_root = aircraft.wing_geom(1,4);
    c_tip  = aircraft.wing_geom(2,4);
    
    % Spanwise arrays for integration
    Y = linspace(0, b2, 50); 
    c_y = interp1([0, b2], [c_root, c_tip], Y);
    
    % 2. THE SHELL MASS (Skin, Ribs, Paint, Trailing Edge)
    % Sandwich composite shell: ~2.8 kg/m^2 per side
    area_density_skin = 2.8; 
    Wing_Area_Half = trapz(Y, c_y); 
    Mass_shell_half = Wing_Area_Half * 2 * area_density_skin; 
    
    % 3. THE STRUCTURAL MASS (Carbon Fiber Spar)
    rho_carbon = 1600; 
    
    % Define the spar geometry along the span
    tc_ratio = 0.1438; 
    t_max_y = c_y .* tc_ratio; 
    % Realistic spar width: 10% of chord
    box_width_y = 0.25 .* c_y; 
    box_height_y = t_max_y .* 0.90;
    
    % Hardcoded spar dimensions
    t_web = 0.003; % 3mm
    t_cap = 0.008; % 8mm
    
    % Calculate cross-sectional area of the spar box at each point
    Area_caps_y = 2 .* (box_width_y .* t_cap);
    Area_webs_y = 2 .* ((box_height_y - 2*t_cap) .* t_web);
    Area_total_spar_y = Area_caps_y + Area_webs_y;
    
    % Integrate Volume along the span and multiply by density
    Volume_spar_half = trapz(Y, Area_total_spar_y);
    Mass_spar_half = Volume_spar_half * rho_carbon;
    
    % 4. Total Wing Weight
    Mass_wing_half = Mass_shell_half + Mass_spar_half;
    W_wing = (Mass_wing_half * 2) * 9.81; % Total wing weight in Newtons
    
    % 5. Non-Lifting Mass (Fuselage, Tail, Pilot)
    M_payload = 90; % Pilot
    M_fuselage_tail = 150; % Realistic empty fuselage/tail mass
    
    % 6. Water Ballast
    W_water = 0; % Set to 0 for dry flight
    
    % 7. Final Total Weight
    W_total = W_wing + (M_payload + M_fuselage_tail)*9.81 + W_water;

end