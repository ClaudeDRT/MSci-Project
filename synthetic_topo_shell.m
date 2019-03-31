function [PosT,surface_type] = synthetic_topo_shell(op_mode,topo_type,sigma_surf,l_surf,H_surf,dx,L_w,L_h,D_off,f_p,L_ang)

%% Generates xyz vectors of synthetic sea ice topography

% Virtual sea ice surface topography with Gaussian, lognormal or fractal
% roughness properties

% Input:
% op_mode: 1 = pulse-limited, 2 = SAR
% topo_type: 1 = Gaussian, 2 = lognormal, 3 = fractal
% sigma_surf = rms roughness height
% l_surf = correlation length
% H_surf = Hurst parameter for fractal surfaces
% dx = grid resolution
% L_w = lead width
% L_h = lead depth
% D_off = lead distance off nadir
% f_p = melt pond fraction

% Output:
% PosT = n x 3 matrix of the xyz surface topography coordinate vertices
% surface_type: 0 = lead/ocean, 1 = sea ice

% (C) Jack Landy, University of Bristol, 2018

%% Create grid

% (Can be modified but these are default grid sizes)
if op_mode < 2
    L = 8000; % across-track diameter of grid, m
    W = 8000; % along-track diameter of grid, m
else
    L = 8000; % across-track diameter of grid, m
    W = 400; % along-track diameter of grid, m
end

[x,y] = meshgrid(-W/2+dx/2:dx:W/2-dx/2,-L/2+dx/2:dx:L/2-dx/2);

%% Generate topography

if topo_type == 1 % Gaussian
    [z,~,~] = rsgene2D_anisotrop(W/dx,L/dx,W,L,sigma_surf,l_surf,l_surf,0);
elseif topo_type == 2 % Lognormal
    [z,~,~] = rsgene2D_anisotrop(W/dx,L/dx,W,L,sigma_surf,l_surf,l_surf,1);
else % Fractal
    [z,~,~] = artificial_surf(sigma_surf,H_surf,W,W/dx,L/dx,(2*pi)/(l_surf*dx));
end

%% Add lead

surface_type = ones(size(z)); % sea ice = 1, lead/ocean = 0

if L_w>0
    if  L_ang ==90
        L_x_min_steps = max([round(L/(2*dx)-(D_off/dx)-L_w/(2*dx)) 1]);
        L_x_max_steps = min([round(L/(2*dx)-(D_off/dx)+L_w/(2*dx))-1 (L/dx)]);
        z(L_x_min_steps:L_x_max_steps,:) = -L_h;
        surface_type(L_x_min_steps:L_x_max_steps,:) = 0;
    elseif L_ang ==0
        L_y_min_steps = max([round(W/(dx*2)-D_off/dx-L_w/(2*dx)) 1]);
        L_y_max_steps = min([round(W/(dx*2)-D_off/dx+L_w/(2*dx))-1 (W/dx)]);
        z(:,L_y_min_steps:L_y_max_steps) = -L_h;
        surface_type(:,L_y_min_steps:L_y_max_steps) = 0;
    elseif L_ang<=45
        a=(round(L_w/(2*dx*sin(L_ang*pi/180))));
        for xi = -a:(L/dx)+a
            line_y = round(tan(L_ang*pi/180)*xi+(W/(2*dx))+(D_off/dx)*cos(L_ang*pi/180)-tan(L_ang*pi/180)*((L/(2*dx))-(D_off/dx)*sin(L_ang*pi/180)));
            line_y_check = round(W/(2*dx) + (D_off/dx)*cos(pi*L_ang/180) + tan(pi*L_ang/180)*(xi - L/(2*dx) + (D_off/dx)*sin(pi*L_ang/180)));
            for yi = 1:(W/dx)
                if yi==line_y
                    L_x_min_steps = max([round(xi-L_w/(2*dx*sin(L_ang*pi/180))) 1]);
                    L_x_max_steps = min([round(xi+L_w/(2*dx*sin(L_ang*pi/180)))-1 (L/dx)]);
                    z(L_x_min_steps:L_x_max_steps, yi) = -L_h;
                    surface_type(L_x_min_steps:L_x_max_steps, yi) =0;
                end
            end
        end
    elseif L_ang>45 & L_ang<90
        b=(round(L_w/(2*dx*cos(pi*L_ang/180))));
        for yi = -b:(W/dx)+b
            line_x = round(L/(2*dx) - (D_off/dx)*sin(L_ang*pi/180) + (1/tan(L_ang*pi/180))*(yi - W/(2*dx) - (D_off/dx)*cos(pi*L_ang/180)));
            for xi = 1:(L/dx)
                if xi==line_x
                    L_y_min_steps = max([round(yi - L_w/(2*dx*cos(pi*L_ang/180))) 1]);
                    L_y_max_steps = min([round(yi + L_w/(2*dx*cos(pi*L_ang/180)))-1 (W/dx)]);
                    z(xi, L_y_min_steps:L_y_max_steps) = -L_h;
                    surface_type(xi, L_y_min_steps:L_y_max_steps) = 0;
                end
            end
        end         
    end
end

%% Add melt ponds

if f_p>0
    [z,surface_type] = add_melt_ponds(z,surface_type,f_p); % topo still referenced to mean height
else
end

%% Finalize

PosT = [x(:) y(:) z(:)];
surface_type = surface_type(:);


end

