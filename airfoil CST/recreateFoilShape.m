% corrected script to test function D_airfoil2 and optimize CST coefficients

clear;
close all;
clc;

global airfoilname
airfoilname = 'e553.dat';

% initial guess coefficients (example)
Au = [1 1 1 1 1];     % upper-surface Bernstein / CST coefficients
Al = [-1 -1 -1 -1 -1];  % lower-surface Bernstein / CST coefficients

% read coordinate file once (format: x y per line)
fid = fopen(airfoilname,'r');
Coor = fscanf(fid, '%g %g', [2 Inf])';
fclose(fid);

% split into X and Y arrays
x_coords = Coor(:,1);
y_coords = Coor(:,2);

% find index of trailing edge (x closest to 1.0).
[~, idx_te] = min(abs(x_coords - 1.0));
if idx_te == 1
    [~, idx_te] = min(abs(x_coords));
end
X_u = x_coords(1:idx_te);
Y_u = y_coords(1:idx_te);
X_l = x_coords(idx_te+1:end);
Y_l = y_coords(idx_te+1:end);

figure;
plot(x_coords,y_coords,'black'); hold on;

% Optimization setup
M = 10;                     % total number of unknowns
x0 = ones(M,1);

options = optimset('Display','Iter');

% run optimization
tic
[x_opt, fval, exitflag, output] = fmincon(@CST_objective, x0, [], [], [], [], [], [], [], options);
toc

% unpack optimized coefficients
Au_opt = x_opt(1:M/2)';
Al_opt = x_opt(M/2+1:end)';

% compute and plot optimized CST
x_sample = linspace(0,1,100)';
[Co_CST_up_opt, ~] = D_airfoil2(Au_opt, Al_opt, x_sample);
[~, Co_CST_low_opt] = D_airfoil2(Au_opt, Al_opt, x_sample);


plot(Co_CST_up_opt(:,1), Co_CST_up_opt(:,2), 'r', 'LineWidth', 1.4);
plot(Co_CST_low_opt(:,1), Co_CST_low_opt(:,2), 'r', 'LineWidth', 1.4);
legend('Airfoil points','CST opt (upper)','CST opt (lower)');
xlabel('x'); ylabel('y');
grid on; ylim([-0.3 0.5]);

% -------------------
% Objective function
% -------------------
function error = CST_objective(x)
    % This function re-reads the airfoil file, splits upper/lower, computes CST
    % surfaces at the same x-locations as the data and returns sum of squared
    % residuals between the y-data and CST y-values.
    global airfoilname
    
    Au = x(1:length(x)/2);
    Al = x(length(x)/2+1:end);
    
    % read the coordinate file
    fid = fopen(airfoilname,'r');
    if fid < 0
        error('Could not open airfoil file: %s', airfoilname);
    end
    Coor = fscanf(fid, '%g %g', [2 Inf])';
    fclose(fid);
    
    x_coords = Coor(:,1);
    y_coords = Coor(:,2);
    
    % find trailing edge index (x closest to 1.0)
    [~, idx_te] = min(abs(x_coords - 1.0));
    
    % split into upper and lower (same logic as main script)
    X_u = x_coords(1:idx_te);
    Y_u = y_coords(1:idx_te);
    X_l = x_coords(idx_te:end);
    Y_l = y_coords(idx_te:end);
    
    % evaluate CST at the data x-locations for upper and lower
    % expecting D_airfoil2 to accept a vector X and return coordinates with same order
    [Co_CST_up, ~] = D_airfoil2(Au, Al, X_u);
    [~, Co_CST_low] = D_airfoil2(Au, Al, X_l);
    
    % If D_airfoil2 returns column vectors where second column is y,
    % compute residuals. If it returns only y-values, adjust accordingly.
    % Here we assume Co_CST_up(:,2) and Co_CST_low(:,2) are y-values.
    
    if size(Co_CST_up,2) < 2
        y_cst_up = Co_CST_up(:);
    else
        y_cst_up = Co_CST_up(:,2);
    end
    if size(Co_CST_low,2) < 2
        y_cst_low = Co_CST_low(:);
    else
        y_cst_low = Co_CST_low(:,2);
    end
    
    % ensure same length (if mismatch, interpolate CST to Y data x-locations)
    if length(y_cst_up) ~= length(Y_u)
        % try to interpolate: assume first column of Co_CST_up is x
        if size(Co_CST_up,2) >= 2
            x_cst_up = Co_CST_up(:,1);
            y_cst_up = interp1(x_cst_up, Co_CST_up(:,2), X_u, 'linear', 'extrap');
        else
            y_cst_up = interp1(linspace(0,1,length(y_cst_up))', y_cst_up, X_u, 'linear', 'extrap');
        end
    end
    if length(y_cst_low) ~= length(Y_l)
        if size(Co_CST_low,2) >= 2
            x_cst_low = Co_CST_low(:,1);
            y_cst_low = interp1(x_cst_low, Co_CST_low(:,2), X_l, 'linear', 'extrap');
        else
            y_cst_low = interp1(linspace(0,1,length(y_cst_low))', y_cst_low, X_l, 'linear', 'extrap');
        end
    end
    
    error_up = Y_u - y_cst_up;
    error_low = Y_l - y_cst_low;
    
    error = sum(error_up.^2) + sum(error_low.^2);
end




%Function D_airfoil2 to transform input Bernstein parameters (shape + class function method) to
%complete Airfoil coordinates

%output [upper surface y coord, lower surface y coord, Class function pos y
%coord, up surf thickness distribution, lw surf thickness distb, camber
%distb] = input (up surf Bernstein parameters, lw surf BS params,
%X-ordinates)
function[Xtu,Xtl,C,Thu,Thl,Cm] = D_airfoil2(Au,Al,X)

x = X(:,1);

N1 = 0.5;   %Class function N1
N2 = 1;     %Class function N2

zeta_u = 0.000;     %upper surface TE gap
zeta_l = -0.000;     %lower surface TE gap


nu = length(Au)-1;
nl = length(Al)-1;

%evaluate required functions for each X-coordinate
for i = 1:length(x)
    
    %calculate Class Function for x(i):
    C(i) = (x(i)^N1)*(1-x(i))^N2;
    
    %calculate Shape Functions for upper and lower surface at x(i)
    Su(i) = 0;  %Shape function initially zero
    for j = 0:nu
        Krnu = factorial(nu)/(factorial(j)*factorial(nu-j));
        Su(i) = Su(i) + Au(j+1)*Krnu*(1-x(i))^(nu-j)*x(i)^(j);
    end
    Sl(i) = 0;  %Shape function initially zero
    for k = 0:nl        
        Krnl = factorial(nl)/(factorial(k)*factorial(nl-k));
        Sl(i) = Sl(i) + Al(k+1)*Krnl*(1-x(i))^(nl-k)*x(i)^(k);
    end
    
    %calculate upper and lower surface ordinates at x(i)
    Yu(i) = C(i)*Su(i) + x(i)*zeta_u;
    Yl(i) = C(i)*Sl(i) + x(i)*zeta_l;
    
    Thu(i) = C(i)*(Su(i)-Sl(i))/2;    %calculate thickness distribution !TE thickness ignored!
    Thl(i) = C(i)*(Sl(i)-Su(i))/2;    %calculate thickness distribution !TE thickness ignored!
    Cm(i) = C(i)*(Su(i)+Sl(i))/2;    %calculate camber distribution !TE thickness ignored!
end

Yust = Yu';
Ylst = Yl';

%assemble airfoil coord matrix
Xtu = [x  Yust];
Xtl = [x  Ylst];

end