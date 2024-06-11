
clear all;
clc;
points =[16, 32, 64, 90, 128, 180];
for n = 1:length(points)
    % Problem parameters
    Lx = 1;           % Length domain in x
    Ly = 1;           % Length domain in y
    Nx = points(4);         % Number of grid points in x
    Ny = points(4);         % Number of grid points in y
    Tfinal = 1;       % Final simulation times
    CFL = 0.2;        % Courant-Friedrichs-Lewy (must be <= 1 for stability)

    % spatial discretization
    dx = Lx / (Nx);
    dy = Ly / (Ny);
    x = linspace(0, Lx, Nx+1);
    y = linspace(0, Ly, Ny+1);
    [X, Y] = meshgrid(x, y);

    % Initial conditions
    k = 10;
    x0 = 0.5;
    y0 = 0.5;
    u0 = exp(-k^2 * ((X-x0).^2 + (Y-y0).^2));

    % Time step
    dt = CFL * dx;
    % Number time steps
    Nt = ceil(Tfinal / dt);
    dt = Tfinal/Nt;

   
    


    % Lax-Wendroff 2-dimension
    u = u0;
    t=0;
    %for n = 1:Nt
    while t < Tfinal
        

        % update time step
        u_half = u;
        for j = 2:Ny

            u_half(1, j) = u(1, j) - 0.5 * CFL * (u(2, j) - u(Nx, j)) + 0.5 * CFL^2 *(u(2, j) - 2*u(1, j) + u(Nx, j));
            u_half(j, 1) = u(j, 1) - 0.5 * CFL * (u(j, 2) - u(j, Nx)) + 0.5 * CFL^2 *(u(j, 2) - 2*u(j, 1) + u(j, Nx));

            for i = 2:Nx
                u_half(i, j) = u(i, j) - 0.5 * CFL * (u(i+1, j) - u(i-1, j)) + 0.5 * CFL^2 * (u(i+1, j) - 2*u(i, j) + u(i-1, j))...
                       - 0.5 * CFL * (u(i, j+1) - u(i, j-1)) + 0.5 * CFL^2 * (u(i, j+1) - 2*u(i, j) + u(i, j-1));
            end


            u_half(Nx+1, j)  = u(Nx+1, j) - 0.5 * CFL * (u(1, j) - u(Nx, j)) + 0.5 * CFL^2 * (u(1, j) - 2*u(Nx+1, j) + u(Nx, j)); 
            u_half(j, Nx+1)  = u(j, Nx+1) - 0.5 * CFL * (u(j, 1) - u(j, Nx)) + 0.5 * CFL^2 * (u(j, 1) - 2*u(j, Nx+1) + u(j, Nx));

        end

        %update boundary conditions
        u_half(1,1) = u(1, 1) - 0.5 * CFL * (u(2, 1) - u(Nx, 1)) + 0.5 * CFL^2 * (u(2, 1) - 2*u(1, 1) + u(Nx, 1)) ...
            - 0.5 * CFL * (u(1, 2) - u(1, Nx)) + 0.5 * CFL^2 * (u(1, Nx) - 2*u(1, 1) + u(1, Nx));

        u_half(1,Nx+1) = u(1, Nx+1) - 0.5 * CFL * (u(1, 1) - u(1, Nx)) + 0.5 * CFL^2 * (u(1, Nx) - 2*u(1, Nx+1) + u(1, 1)) ...
            - 0.5 * CFL * (u(2, Nx+1) - u(Nx, Nx+1)) + 0.5 * CFL^2 * (u(Nx, Nx+1) - 2*u(1, Nx+1) + u(2, Nx+1));

        u_half(Nx+1, 1) = u(Nx+1, 1) - 0.5 * CFL * (u(1, 1) - u(Nx, 1)) + 0.5 * CFL^2 * (u(1, 1) - 2*u(Nx+1, 1) + u(Nx, 1)) ...
            - 0.5 * CFL * (u(Nx+1,2) - u(Nx+1, Nx)) + 0.5 * CFL^2 * (u(Nx+1, Nx) - 2*u(Nx+1, Nx) + u(Nx+1, 2));

        u_half(Nx+1, Nx+1) = u(Nx+1, Nx+1) - 0.5 * CFL * (u(1, Nx+1) - u(Nx, Nx+1)) + 0.5 * CFL^2 * (u(1, Nx+1) - 2*u(Nx+1, Nx+1) + u(Nx, Nx+1)) ...
            - 0.5 * CFL * (u(Nx+1, 1) - u(Nx+1, Nx)) + 0.5 * CFL^2 * (u(Nx+1, 1) - 2*u(Nx+1, Nx+1) + u(Nx+1, Nx));    

        % Actualizar la solución en el tiempo

        u = u_half;
        % Update time
        t = t + dt;



        % Plot the solution at each time step (optional)
%         surf(X, Y, u);
%         axis([0 Lx 0 Ly -1 1]);
%         title(['Time step: ' num2str(n) '/' num2str(Nt)]);
%         xlabel('x');
%         ylabel('y');
%         zlabel('u');
%         drawnow;
    end
                 Error = u-u0;
        numErrorsL2(n)  = sqrt(sum(sum((Error).^2)) / (Nx^2));
        surf(X, Y, u);
        axis([0 Lx 0 Ly -1 2]);
        xlabel('x');
        ylabel('y');
        zlabel('u');
        ax = gca;
        ax.FontSize = 15;  % Adjust the number to your preferred size
        drawnow;

end




