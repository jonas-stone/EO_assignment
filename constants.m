function c = constants()

%% CONSTANTS
% function for the physical constants that will not change

% Inputs: 
% Outputs: c
c = struct();

%% Flight Conditions =======================================================
c.altitude = 1500;                       % Flying altitude              [m]

%% Physical constants ======================================================
c.g                       = 9.81;        % Gravity                  [m/s^2]
c.R                       = 287;         % Air gas constant      [J/(kg K)]
c.gamma                   = 1.4;         % Air ratio of specific heats  [/]

%% Atmospheric constants ===================================================
c.T0                      = 288.15;      % Sea-level temperature        [K]
c.P0                      = 101325;      % Sea-level pressure          [Pa]
c.h_trop                  = 11000;       % Tropopause altitude          [m]
c.T_trop                  = 216.65;      % Tropopause temperature       [K]
c.P_trop                  = 22632;       % Pressure  at 11 km          [Pa]
c.lapse                   = 0.0065;      % Temperature lapse rate     [K/m]

%% Thermodynamic constants =================================================
c.mu0                     = 1.716e-5;    % Ref. dynamic viscosity    [Pa s]
c.Tref                    = 273.15;      % Ref. temperature             [K]
c.S                       = 110.4;       % Sutherland constant          [K]

%% Material properties =====================================================
c.material.E    = 6969e6;     % Young's modulus             [Pa]
c.material.rho  = 6969;        % Density                 [kg/m^3]
c.material.Ft_y = 6969;       % Tensile yield strength      [Pa]
c.material.Fc_y = 6969;       % Compressive yield strength  [Pa]


end