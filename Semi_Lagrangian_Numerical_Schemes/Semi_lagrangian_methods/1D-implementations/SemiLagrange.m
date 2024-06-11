clear all;
clc;

% Semi-Lagrangian Scheme for 1D Advection Equation

points =[64,90,128,180,256,360,512,720];
count=1;
for i = 1:1
    % Parameters
    L = 1;           % Length of the domain
    Nx = points(5);  % Number of grid points
    dx = L / (Nx);   % Grid spacing
    c = 1;           % Advection speed
   
    
    % Time-stepping loop
    tFinal = 15.0;  % Final simulation time
    t = 0;
    
    CFL = 0.8;  % Desired CFL number 
     dt = CFL * dx ;
    Nt = ceil(L / dt);
    dt = L / Nt;

     Nt = floor(tFinal/dt);
     dt = tFinal/Nt;
    
    % Initial condition
    k = 25;
    
    x0 = 0.5;
    u0 =@(x) exp(-k^2 * (x-x0).^2);
    
    % Initialize grid and solution
    x = linspace(0, L, Nx);
    u = u0(x);
    
    
    for j = 1:Nt
        % Semi-Lagrangian update
         x_back = mod(x - c*dt, L);  % Backward particle tracing
        
        %--- SPLINE/LINEAR INTERPOLATION ------
%       u1 = interp1(x, u, x_back, 'linear', 0);  % Spline interpolation
        
        %--- LAGRANGE INTERPOLATION ------
       [idx, ~] = knnsearch(x',x_back','K',5);
        for j=1:Nx
            u1(j) = lagrange_interpolation(x(idx(j,:)), u(idx(j,:)), x_back(j));
        end    
        
        u = u1;
        
        % Plot the solution
%         plot(x, u, '-o', 'LineWidth', 2);
%         title(['Time = ', num2str(t)]);
%         xlabel('x');
%         ylabel('u');
%         axis([0, L, 0, 1]);
%         drawnow;
%         
        % Update time
        t = t + dt;
        
    end
 
end





% Mostrar el resultado final
figure;
plot(x, u0(x), '-o', x, u, '-o');
axis([0 L -0.1 1.1])
legend('Initial condition', 'Final solution');
xlabel('x');
ylabel('u(x,T)');
ax = gca;
ax.FontSize = 15;  % Adjust the number to your preferred size



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


