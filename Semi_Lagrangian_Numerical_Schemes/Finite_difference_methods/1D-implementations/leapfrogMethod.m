

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
    
    % Leapfrog method
    u = u0;
    u_half = zeros(size(u));
    utp1 = zeros(size(u));
    
    % Calcular el paso en el tiempo inicial usando un paso de Euler
    for i = 2:Nx
        u_half(i) = u(i) -  CFL * (u(i) - u(i-1));
    end
    
    % Condiciones de contorno periódicas
    u_half(1) = u(1) -  CFL * (u(1) - u(Nx));
    u_half(Nx+1) = u_half(1);
    
    for n = 2:Nt
        % Calculate time step 
        for i = 2:Nx
            utp1(i) = u(i) - CFL * (u_half(i+1) - u_half(i-1));
        end
    
        % Boundary conditions
        utp1(1) = u(1) - CFL * (u_half(2) - u_half(Nx));
        utp1(Nx+1) = utp1(1);
    
        %update solution
        u = u_half;
        u_half = utp1;
        
    
    
        % Plot solution (optional)
%         plot(x, utp1, '-o');
%         axis([0 L -1 1]);
%         title(['Paso de tiempo: ' num2str(n) '/' num2str(Nt)]);
%         drawnow;
    end

end


figure;
plot(x, u0, '-o', x, utp1, '-o');
axis([0 L -0.5 1.3])
legend('Initial condition', 'Final solution');xlabel('x');
ylabel('u(x,T)');
ax = gca;
ax.FontSize = 15;  % Adjust the number to your preferred size

