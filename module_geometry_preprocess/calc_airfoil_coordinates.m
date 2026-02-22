function airfoil_coords = calc_airfoil_coordinates(Au, Al)
    %{
    %% CALC AIRFOIL COORDINATES FUNCTION

    %% Inputs:
    - Au
    - Al
    %% Outputs:
    - airfoil_coords
    
    %% Notes:
    Turning CST coefficients into a .dat file with the airfoil coordinates.    
    %}

    % Chord points sample
    x_sample = linspace(0, 1, 1000)';

    % Get upper and lower CST coordinates
    [xz_CST_up, xz_CST_lo] = D_airfoil2(Au, Al, x_sample);

    % Create airfoil surface points: flipped upper surface, then lower 
    % surface
    airfoil_coords = [flipud(xz_CST_up); xz_CST_lo];

end