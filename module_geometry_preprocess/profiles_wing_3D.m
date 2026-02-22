function wing_profiles = profiles_wing_3D(Au_initial, Al_initial, wing_geom, n_divs, plot_on)
	%{
	UPDATED TO JONAS
AI HAS BEEN USED TO HELP WITH PLOTTING
	%}
chords = wing_geom(:,4);
xyz_le = wing_geom(:,1:3);
twists = wing_geom(:,5);

wing_profiles = struct();
wing_profiles.in.up = ones(n_divs,3);
wing_profiles.in.lo = ones(n_divs,3);
% wing_profiles.mid.up = ones(n_divs,3);
% wing_profiles.mid.lo = ones(n_divs,3);
wing_profiles.out.up = ones(n_divs,3);
wing_profiles.out.lo = ones(n_divs,3);

for i = 1:length(chords)
    x_sample = linspace(0,1,n_divs)';
    [xz_upper, xz_lower] = D_airfoil2(Au_initial, Al_initial, x_sample);

    xz_upper = xz_upper * chords(i);
    xz_lower = xz_lower * chords(i);

    theta = deg2rad(twists(i));
    R = [cos(theta), -sin(theta); 
         sin(theta),  cos(theta)];

    xz_upper = (R * xz_upper')';
    xz_lower = (R * xz_lower')';

    if i == 1
        wing_profiles.in.up = [xz_upper(:,1) + xyz_le(i,1), ones(n_divs,1)*xyz_le(i,2), xz_upper(:,2) + xyz_le(i,3)];
        wing_profiles.in.lo = [xz_lower(:,1) + xyz_le(i,1), ones(n_divs,1)*xyz_le(i,2), xz_lower(:,2) + xyz_le(i,3)];
    elseif i == length(chords)
        wing_profiles.out.up = [xz_upper(:,1) + xyz_le(i,1), ones(n_divs,1)*xyz_le(i,2), xz_upper(:,2) + xyz_le(i,3)];
        wing_profiles.out.lo = [xz_lower(:,1) + xyz_le(i,1), ones(n_divs,1)*xyz_le(i,2), xz_lower(:,2) + xyz_le(i,3)];
    else
        % wing_profiles.mid.up = [xz_upper(:,1) + xyz_le(i,1), ones(n_divs,1)*xyz_le(i,2), xz_upper(:,2) + xyz_le(i,3)];
        % wing_profiles.mid.lo = [xz_lower(:,1) + xyz_le(i,1), ones(n_divs,1)*xyz_le(i,2), xz_lower(:,2) + xyz_le(i,3)];
    end
    
end


%% Visualize wing in 3D 
% Number of spanwise panels for smooth surface
n_span = 80;

% Extract coordinates for upper and lower surfaces (n_divs x 3)
XU = [wing_profiles.in.up(:,1),  wing_profiles.out.up(:,1)];
YU = [wing_profiles.in.up(:,2),  wing_profiles.out.up(:,2)];
ZU = [wing_profiles.in.up(:,3),  wing_profiles.out.up(:,3)];

XL = [wing_profiles.in.lo(:,1),  wing_profiles.out.lo(:,1)];
YL = [wing_profiles.in.lo(:,2),  wing_profiles.out.lo(:,2)];
ZL = [wing_profiles.in.lo(:,3),  wing_profiles.out.lo(:,3)];

% Interpolate spanwise to create smooth surface
xi = linspace(1,2,n_span);
XU_s = interp1(1:2, XU', xi, 'linear')';
YU_s = interp1(1:2, YU', xi, 'linear')';
ZU_s = interp1(1:2, ZU', xi, 'linear')';

XL_s = interp1(1:2, XL', xi, 'linear')';
YL_s = interp1(1:2, YL', xi, 'linear')';
ZL_s = interp1(1:2, ZL', xi, 'linear')';

if plot_on == 1
    % Plot
    figure('Color','w'); hold on; grid on;
    
    % Upper surface
    surf(XU_s, YU_s, ZU_s, 'FaceColor',[0.6 0.8 1], 'EdgeColor','none', 'FaceAlpha',0.2);
    
    % Lower surface
    surf(XL_s, YL_s, ZL_s, 'FaceColor',[1 0.9 0.8], 'EdgeColor','none', 'FaceAlpha',0.2);
    
    % Close root and tip by filling between upper and lower cross-sections
    % Root (first column)
    root_x = [XU_s(:,1); flip(XL_s(:,1))];
    root_y = [YU_s(:,1); flip(YL_s(:,1))];
    root_z = [ZU_s(:,1); flip(ZL_s(:,1))];
    fill3(root_x, root_y, root_z, [0.7 0.85 0.95], 'EdgeColor','none');
    
    % Tip (last column)
    tip_x = [XU_s(:,end); flip(XL_s(:,end))];
    tip_y = [YU_s(:,end); flip(YL_s(:,end))];
    tip_z = [ZU_s(:,end); flip(ZL_s(:,end))];
    fill3(tip_x, tip_y, tip_z, [0.95 0.9 0.8], 'EdgeColor','none');
    
    % Leading and trailing edge lines (use mean of upper/lower for trailing)
    LE_x = 0.5*(XU_s(1,:) + XL_s(1,:));
    LE_y = 0.5*(YU_s(1,:) + YL_s(1,:));
    LE_z = 0.5*(ZU_s(1,:) + ZL_s(1,:));
    plot3(LE_x, LE_y, LE_z, 'k-', 'LineWidth',1.5);
    
    TE_x = 0.5*(XU_s(end,:) + XL_s(end,:));
    TE_y = 0.5*(YU_s(end,:) + YL_s(end,:));
    TE_z = 0.5*(ZU_s(end,:) + ZL_s(end,:));
    plot3(TE_x, TE_y, TE_z, 'k-', 'LineWidth',1.5);
    
    % Optional: also plot section outlines (three original sections)
    cols = {'b','r'};
    sections = {'in','out'};
    for k = 1:2
        up = wing_profiles.(sections{k}).up;
        lo = wing_profiles.(sections{k}).lo;
        plot3(up(:,1), up(:,2), up(:,3), '-', 'Color', cols{k}, 'LineWidth', 1);
        plot3(lo(:,1), lo(:,2), lo(:,3), '--', 'Color', cols{k}, 'LineWidth', 1);
    end
    
    % Aesthetics
    axis equal;
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title('3D Wing Visualization (upper & lower surfaces)');
    lighting gouraud;
    camlight headlight;
    material dull;
    view(3);
    rotate3d on;
end
end