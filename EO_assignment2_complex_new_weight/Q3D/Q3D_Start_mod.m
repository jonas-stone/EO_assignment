function [Res, LD, L] = Q3D_Start_mod(aircraft,aero)

cd("Q3D\");

% Wing planform geometry 
%                x    y     z   chord(m)    twist angle (deg) 
AC.Wing.Geom = aircraft.wing_geom;

% Wing incidence angle (degree)
AC.Wing.inc  = aircraft.wing_geom(1,4);

% Airfoil coefficients input matrix
AC.Wing.Airfoils   = aircraft.airfoils;
                  
AC.Wing.eta = [0;1];  % Spanwise location of the airfoil sections

% Viscous vs inviscid
AC.Visc  = 1;              % 0 for inviscid and 1 for viscous analysis
AC.Aero.MaxIterIndex = 60;    %Maximum number of Iteration for the
                                %convergence of viscous calculation

% Flight Condition
AC.Aero.V     = aero.V;         % flight speed (m/s)
AC.Aero.rho   = aero.rho;       % air density  (kg/m3)
AC.Aero.alt   = aero.alt;       % flight altitude (m)
AC.Aero.Re    = aero.Re;        % reynolds number (bqased on mean aerodynamic chord)
AC.Aero.M     = aero.M;         % flight Mach number 
% AC.Aero.CL    = aero.CL;     % angle of attack -  comment this line to run the code for given cl 
AC.Aero.Alpha = aero.Alpha;     % angle of attack -  comment this line to run the code for given cl 

%% 
tic
Res = Q3D_solver(AC);
toc

if ~isnan(Res.CDwing) 
    % --- NEW: ADD FUSELAGE & TAIL DRAG ---
    % Equivalent parasite area of a modern sailplane fuselage & tail (m^2)
    f_area = 0.09; 
    
    % Convert parasite area to a Drag Coefficient based on the current wing area
    CD_fuse_tail = f_area / aircraft.S; 
    
    % Sum the wing drag and the non-lifting drag
    Total_CD = Res.CDwing + CD_fuse_tail;
    
    % Calculate true Aircraft L/D
    LD = Res.CLwing / Total_CD;
    
    % Lift force calculation (remains unchanged, using CLwing)
    L = Res.CLwing * 1/2 * aero.V^2 * aero.rho * aircraft.S;
    % -------------------------------------
else
    LD = NaN;
    L = NaN;
end

cd("..\");

end