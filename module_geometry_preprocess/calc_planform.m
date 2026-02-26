function aircraft = calc_planform(b2,c_root,c_tip,twist_tip)
	%{
	UPDATED TO JONAS
    %}
aircraft.b2 = b2;
aircraft.c_root = c_root;
aircraft.c_tip = c_tip;
aircraft.twist_tip = twist_tip;
%% -------------------------------
%  Constant Angles Definitons (deg)
aircraft.LE_sweep   = 5;
aircraft.dihedral   = 5;
aircraft.root_incidence = 0;

%% -------------------------------
%  Spanwise stations
y_root = 0;
y_tip  = aircraft.b2;

%% -------------------------------
%  Leading-edge x-positions
x_le_root = 0;
x_le_tip  =  y_tip * tand(aircraft.LE_sweep);

%% -------------------------------
%  Dihedral (z positions)
% -------------------------------
z_root = 0;
z_tip  = y_tip * tand(aircraft.dihedral);

%% -------------------------------
%  Assemble wing_geom matrix
% Format: [x_le, y_le, z_le, chord, twist]
% -------------------------------
aircraft.wing_geom = [
    x_le_root, y_root, z_root, c_root, aircraft.root_incidence;  % Root
    x_le_tip,  y_tip,  z_tip,  c_tip,  aircraft.twist_tip    % Tip
];

%% -------------------------------
%  8. Wing reference area (half + full)
aircraft.S_half = (c_root + c_tip) / 2 * aircraft.b2; % A = (c_root + c_tip)/2 * span
aircraft.S      = 2 * aircraft.S_half;

%% -------------------------------
%  9. Full span and Mean Aerodynamic Chord
aircraft.b = 2 * aircraft.b2;
% Taper ratio
aircraft.lambda = c_tip / c_root;

aircraft.AR = (aircraft.b^2) / aircraft.S;

% MAC length
aircraft.MAC = (2/3) * c_root * ((1 + aircraft.lambda + aircraft.lambda^2) / (1 + aircraft.lambda));

%% -------------------------------
%  10. Airfoil Definition
%                      | ->     upper curve coeff.                <-|   | ->       lower curve coeff.       <-| 
aircraft.airfoils   = [0.2171    0.3450    0.2975    0.2685    0.2893  -0.1299   -0.2388   -0.1635   -0.0476    0.0797;
                       0.2171    0.3450    0.2975    0.2685    0.2893  -0.1299   -0.2388   -0.1635   -0.0476    0.0797];

end
