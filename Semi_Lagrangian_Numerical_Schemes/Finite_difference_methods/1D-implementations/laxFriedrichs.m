

points =[64,90,128,180,256,360,512,720];
count=1;
for i = 1:1
 
    % Problem parameters
    L = 1;           % Length domain
    Nx = points(5);        % Number of grid points
    Tfinal = 5.0;    % Final simulation time
    CFL = 0.8;       % Courant-Friedrichs-Lewy (must be <= 1 for stability)
    
    % Spatial discretications
    dx = L / Nx;
    x = linspace(0, L, Nx+1);
    
    % Initial conditions
    k=10;
    
    x0 = 0.5;
    u0 = exp(-k^2 * (x-x0).^2);
    
    % Time-step
    dt =  CFL * dx;
    
    % Number of time steps
    Nt = ceil(Tfinal / dt);
    
    % Lax-Firedrichs methods
    u = u0;
    for n = 1:Nt
        % Calculate new time step
        u_half = zeros(size(u));
        for i = 2:Nx
            u_half(i) = 0.5 * (u(i+1) + u(i-1)) - 0.5 * CFL * (u(i+1) - u(i-1));
        end
       
        % Periodic boundary conditions
        u_half(1) = 0.5 * (u(2) + u(Nx)) - 0.5 * CFL * (u(2) - u(Nx));
        u_half(Nx+1) = 0.5 * (u(1) + u(Nx)) - CFL/2 *(u(1)-u(Nx));
       
        % Update solution
        u = u_half;
       
    %     %Plot solution (optional)
    %     plot(x, u, '-o');
    %     axis([0 L -1 1]);
    %     title(['Paso de tiempo: ' num2str(n) '/' num2str(Nt)]);
    %     drawnow;
    end
    


end



% Mostrar el resultado final
figure;
plot(x, u0, '-o', x, u, '-o');
axis([0 L -0.5 1.3])
legend('Initial condition', 'Final solution');
xlabel('x');
ylabel('u(x,T)');
ax = gca;
ax.FontSize = 15;  % Adjust the number to your preferred size


