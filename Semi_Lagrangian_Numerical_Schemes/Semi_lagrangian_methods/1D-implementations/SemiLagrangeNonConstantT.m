% Semi-Lagrangian Scheme for 1D Advection Equation
clear all;
clc;

points = [64, 90, 128, 180, 256, 360, 512, 720];
for i = 1:1
    % Parameters
    L = pi; %1       % Length of the domain
    Nx = points(6);  % Number of grid points
    dx = L / Nx;     % Grid spacing
    
    % Time-stepping loop
    tFinal = 8.0;   % Final simulation time
    t = 0;
    tdata = t;
    
    CFL = 0.8;       % Desired CFL number 
    dt = CFL * dx;
    Nt = floor(tFinal/dt);
    dt = tFinal/Nt;
    
    % Initial condition
    k = 25;
    x0 = 0.5;
    u0 = @(x) exp(-144*(x-1.5).^2);
    
    Vt = @(t)  pi/2 * sin(pi/2 * t + pi/2);
    % Initialize grid and solution
    x = linspace(0, L, Nx);
    u = u0(x);
    
    tplot = 0.05;
    plotgap = round(tplot/dt); 
    dt = tplot/plotgap;
    nplots = round(tFinal/tplot);
    data = [u0(x); zeros(nplots, Nx)];
    
    for z = 1:nplots
        for n=1:plotgap
            % Calculate velocity at current time
            %v_t = pi/2 * sin(pi/2 * t + pi/2);

            % RK4 for time integration of the velocity field
            k1 = dt * Vt(t);
            k2 = dt * Vt(t + 0.5 * dt);
            k3 = dt * Vt(t + 0.5 * dt);
            k4 = dt * Vt(t + dt);

            % Semi-Lagrangian update
            x_back = mod(x - (k1 + 2*k2 + 2*k3 + k4)/6, L);  % Backward particle tracing

%             u1 = interp1(x, u, x_back, 'spline', 0);  % Spline interpolation
            
             [idx, ~] = knnsearch(x',x_back','K',3);
             for j=1:Nx
                 u1(j) = lagrange_interpolation(x(idx(j,:)), u(idx(j,:)), x_back(j));
             end   

            u = u1;

            % Update time
            t = t + dt;
            
          % Plot the solution
%           plot(x, u, '-o', 'LineWidth', 2);
%             title(['Time = ', num2str(t)]);
%             xlabel('x');
%             ylabel('u');
%             axis([0, L, 0, 1]);
%             drawnow;
        end
         data(z+1,:) = u; tdata = [tdata;t];

    end
    

end

% Plot the final result

  waterfall(x,tdata,data);
  view(10,70);
  colormap([0 0 0]);
    axis([0 pi 0 tFinal 0 5]);
    ylabel t;
    zlabel u;
    xlabel x;
    grid off;



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

