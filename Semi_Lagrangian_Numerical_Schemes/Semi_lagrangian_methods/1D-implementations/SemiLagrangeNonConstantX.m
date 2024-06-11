% Semi-Lagrangian Scheme for 1D Advection Equation
clear all;
clc;

points = [64, 90, 128, 180, 256, 360, 512, 720];
count = 1;
for i = 1:length(points)
    % Parameters
%     L = pi; %1           % Length of the domain
    Nx = points(i);  % Number of grid points
    h = 2*pi /Nx;     % Grid spacing
    x = h*(1:Nx);
    dt = h/4;
    
    % Time-stepping loop
    tFinal = 8.0;   % Final simulation time
    t = 0;
    tplot = 0.15;
    plotgap = round(tplot/dt); 
    dt = tplot/plotgap;
    nplots = round(tFinal/tplot);
    tdata = t;
    
    % Initial condition
    k = 25;
    x0 = 0.5;
    u0 = @(x) exp(-100*(x-1).^2);
    data = [u0(x); zeros(nplots, Nx)];
    
    V = @(x)  0.2 + sin(x-1).^2;
    % Initialize grid and solution
    u = u0(x);
    
    for j = 1:nplots
        % Calculate velocity at current time
        %v_t = pi/2 * sin(pi/2 * t + pi/2);
        for n = 1:plotgap
        
            % RK4 for time integration of the velocity field
            k1 = dt * V(x);
            k2 = dt * V(x + 0.5 * k1);
            k3 = dt * V(x + 0.5 * k2);
            k4 = dt * V(x + k3);

            % Semi-Lagrangian update
            x_back = mod(x - (k1 + 2*k2 + 2*k3 + k4)/6, x(end));  % Backward particle tracing

            %----- SPLINE/LINEAR INTERPOLATION ------
             u1 = interp1(x, u, x_back, 'spline', 0);  % Spline interpolation
             
             
            %----- LAGRANGE INTERPOLATION ------
%              [idx, ~] = knnsearch(x',x_back','K',5);
%              for z=1:Nx
%                  u1(z) = lagrange_interpolation(x(idx(z,:)), u(idx(z,:)), x_back(z));
%              end   

            u = u1;

            % Update time
            t = t + dt;
         
          %plot(x, u, '-o', 'LineWidth', 2);
%         title(['Time = ', num2str(t)]);
%         xlabel('x');
%         ylabel('u');
%         axis([0, L, 0, 1]);
%         drawnow;
        end
        data(j+1,:) = u; tdata = [tdata;t];

    end
    

    
    waterfall(x,tdata,data);
    view(10,70);
    colormap([0 0 0]);
    axis([0 2*pi 0 tFinal 0 5]);
    ylabel t;
    zlabel u;
    xlabel x;
    grid off;
    
    % Increment counter
    count = count + 1;
end





function interpolated_values = lagrange_interpolation(x, y, x_interpolate)
    % x: x-coordinates of data points
    % y: y-coordinates of data points
    % x_interpolate: x-values where interpolation is desired
    
    n = length(x);
    interpolated_values = zeros(size(x_interpolate));

    for i = 1:length(x_interpolate)
        % Lagrange basis polynomial
        L = ones(size(x));
        for j = 1:n
            for k = 1:n
                if k ~= j
                    L(j) = L(j) * (x_interpolate(i) - x(k)) / (x(j) - x(k));
                end
            end
        end

        % Interpolated value at x_interpolate(i)
        interpolated_values(i) = sum(y .* L);
    end
end

