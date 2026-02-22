function aircraft = calc_planform(b2,c_root,c_tip,twist_tip)
	%{
	UPDATED TO JONAS
    %}

%% -------------------------------
%  1. Constant Angles Definitons (deg)
LE_sweep   = 5;
dihedral   = 5;
root_incidence = 0;

%% -------------------------------
%  2. Spanwise stations
y_root = 0;
y_tip  = b2;

%% -------------------------------
%  3. Leading-edge x-positions
x_le_root = 0;
x_le_tip  =  y_tip * tand(LE_sweep);

%% -------------------------------
%  4. Trailing-edge x-positions
% Root TE
x_te_root = x_le_root + c_root;

% Outboard TE implied by tip chord
x_te_tip = x_le_tip + c_tip*cosd(twist_tip);

%% -------------------------------
%  6. Dihedral (z positions)
% -------------------------------
z_root = 0;
z_tip  = y_tip * tan(dihedral);

%% -------------------------------
%  7. Assemble wing_geom matrix
% Format: [x_le, y_le, z_le, chord, twist]
% -------------------------------
aircraft.wing_geom = [
    x_le_root, y_root, z_root, c_root, root_incidence;  % Root
    x_le_tip,  y_tip,  z_tip,  c_tip,  twist_tip    % Tip
];

%% -------------------------------
%  8. Wing reference area (half + full)
aircraft.S_ref_half = (c_root + c_tip) / 2 * b2; % A = (c_root + c_tip)/2 * span
aircraft.S_ref      = 2 * aircraft.S_ref_half;

%% -------------------------------
%  9. Full span and Mean Aerodynamic Chord
aircraft.b = 2 * aircraft.b2;
% Taper ratio
lambda = c_tip / c_root;
% MAC length
aircraft.MAC = (2/3) * c_root * ((1 + lambda + lambda^2) / (1 + lambda));

%% -------------------------------
%  10. Airfoil Definition
%                      | ->     upper curve coeff.                <-|   | ->       lower curve coeff.       <-| 
aircraft.airfoils   = [0.2171    0.3450    0.2975    0.2685    0.2893  -0.1299   -0.2388   -0.1635   -0.0476    0.0797;
                       0.2171    0.3450    0.2975    0.2685    0.2893  -0.1299   -0.2388   -0.1635   -0.0476    0.0797];

end
