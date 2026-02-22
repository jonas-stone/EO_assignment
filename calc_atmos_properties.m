function [a, rho, v, mach, mu] = calc_atmos_properties(h, speed_val, inputType)
    %{
    %% CALC ATMOS PROEPRTIES FUNCTION
    
    %% Inputs: 
    -h
    -speed_val
	-inputType

    %% Outputs: 
    - a
    - rho
	- v
    - mach
    - mu

    %% Notes: 
    Calculates atmospheric properties and and other props for a h and v
    %}
    %% Extracting constants
    if nargin < 3
        inputType = 'v';
    end
    c = constants();

    % giving the extracted constants easier names to code with
    g     = c.g;          % Gravity                       [m/s^2]
    R     = c.R;          % Specific gas constant for air [J/(kg·K)]
    gamma = c.gamma;      % Adiabatic index               [-]
    T0     = c.T0;        % Sea-level temperature         [K]
    P0     = c.P0;        % Sea-level pressure            [Pa]
    h_trop = c.h_trop;    % Tropopause altitude           [m]
    T_trop = c.T_trop;    % Tropopause temperature        [K]
    P_trop = c.P_trop;    % Pressure at tropopause        [Pa]
    L      = c.lapse;     % Lapse rate                    [K/m]
    mu0  = c.mu0;         % Reference viscosity           [Pa·s]
    Tref = c.Tref;        % Reference temperature         [K]
    S    = c.S;           % Sutherland constant           [K]

    %% Temperature and Pressure calculations
    if h <= h_trop
        T = T0 - L * h;
        exponent = g / (R * L);
        P = P0 * (T / T0)^exponent;
    else
        T = T_trop;
        delta_h = h - h_trop;
        P = P_trop * exp(-g * delta_h / (R * T_trop));
    end
    
    %% Other properties 
    % Speed of sound
    a = sqrt(gamma * R * T);

    % Density
    rho = P / (R * T);

    % % True airspeed from Mach number
    % mach = v / a;

    % Dynamic viscosity (Sutherland's law)
    mu = mu0 * (T / Tref)^(3/2) * (Tref + S) / (T + S);

    switch lower(inputType)
        case 'mach'
            mach = speed_val;
            v    = mach .* a;  % Calculate TAS
        case {'v'}
            v    = speed_val;
            mach = v ./ a;     % Calculate Mach
        otherwise
            error('Invalid inputType. Use ''Mach'' or ''v''.');
    end	   
end

