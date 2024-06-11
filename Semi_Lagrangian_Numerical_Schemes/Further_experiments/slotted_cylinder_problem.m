% -- ISKE EXPERIMENT 2

clear all;
clc;
% Parameters
Lx = 1;           % Length of the domain in x-direction
Ly = 1;           % Length of the domain in y-direction
Nx =  90;         % Number of grid points in x-direction %40
Ny =  90;         % Number of grid points in y-direction %40 
dx = Lx / Nx;     % Grid spacing in x-direction
dy = Ly / Ny;     % Grid spacing in y-directio

% Time parameters
% tFinal = 6.28;     % Final simulation time %6.28 with constant velocity
tFinal = 5.36;     % Final simulation time %5.36 with aceleration
CFL = 0.5;        % Desired CFL number

% Initialize grid and solution
x = linspace(-0.5, 0.5, Nx);
y = linspace(-0.5, 0.5, Ny);

[X,Y]= meshgrid(x,y);
Z = slotted_cylinder(X,Y);
 Z0 = Z;
 
Nt = 800;
dt = tFinal/Nt;



t = 0;


for i= 1:Nt
    
    
    % Calculate velocity components at each grid point
        [Vx, Vy] = calculate_velocity(X, Y);
    %   [Vx, Vy] = calculate_velocity_constant(X, Y);


    % Semi-Lagrangian update using 
    k1x = dt * Vx;
    k1y = dt * Vy;


    x_back = X - (k1x);  % Backward particle tracing
    y_back = Y - (k1y);  % Backward particle tracing
    

   
    
%     %----------- RBF INTERPOLATION --------
%      %% este funciona pero resultado regular, parecido a interp2 spline
%              pd=4;
%              ns=2*nchoosek(pd+2,2);
%              [ind,~] = knnsearch([X(:) Y(:)],[x_back(:) y_back(:)],'K',ns);
%              for j =1:numel(X)
%                  W = InterpRBF([X(ind(j,:))' Y(ind(j,:))'],[x_back(j) y_back(j)],3,pd);
%                  Z1(j) = sum(W.*Z(ind(j,:))'); % approx. using RBFs
%              end
% %     
%              Z = reshape(Z1,[Nx Nx]);
    

    
    
    %-----------------SPLINE INTERPOLATION -----------------
    
    Z = interp2(X, Y, Z, x_back, y_back, 'spline', 0);  % 2D interpolation
     

    
    
    % Update time
    t = t + dt;
    
%         % Plot the solution at each time step (optional)
%         surf(X, Y, Z);
%         axis([-0.5 0.5 -0.5 0.5 -0.1 1.3]);
%         title(['Time = ', num2str(t)]);
%         xlabel('x');
%         ylabel('y');
%         zlabel('u');
%         drawnow;
%         

    
end

    figure;
         surf(X, Y, Z);
        axis([-0.5 0.5 -0.5 0.5 -0.5 1.5]);
        title(['Time = ', num2str(t)]);
        xlabel('x');
        ylabel('y');
        zlabel('u');
        drawnow;
        
        figure;
        
          surf(X, Y, Z0);
        axis([-0.5 0.5 -0.5 0.5 -0.5 1.5]);
        title(['Time = ', num2str(t)]);
        xlabel('x');
        ylabel('y');
        zlabel('u');
        drawnow;


    


function Z = slotted_cylinder(X,Y)

    X0=0; Y0=0.25;
    R = 0.15;
    d = 0.06;
    L = 0.22;

    Z = zeros(size(X));
    Z((X - X0).^2 + (Y - Y0).^2 < R^2) = 1;

    Z(find((Y > Y0 - d/2) .* (Y < Y0 + d/2) .* (X > X0) .* (X < X0 + L))) = 0;
end




function [vx, vy] = calculate_velocity(x, y)
    % Define azimuth angle phi
    %     phi = atan2(y, x);
    phi(x >= 0) = atan2(-y(x >= 0), x(x >= 0));
    phi(x < 0) = atan2(x(x < 0), y(x < 0)) + pi/2;

    % Define velocity components based on conditions
    a = zeros(size(x));
    a(y <= 0) = 0.5 * sin(2 * phi(y <= 0) - pi/2) + 3/2;
    a(y >= 0) = 1;

    % Calculate velocity components vx and vy
    vx = y .* a;
    vy = -x .* a;
end

