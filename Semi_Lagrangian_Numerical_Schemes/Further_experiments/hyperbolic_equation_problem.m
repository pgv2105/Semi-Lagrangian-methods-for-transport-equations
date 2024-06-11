clear all;
clc;

% Parameters
Lx = 1;           % Length of the domain in x-direction
Ly = 1;           % Length of the domain in y-direction
Nx =  128;         % Number of grid points in x-direction
Ny =  128;         % Number of grid points in y-direction
dx = 1 / Nx;     % Grid spacing in x-direction
dy = 1 / Ny;     % Grid spacing in y-direction

% Time parameters
tFinal = 1;       % Final simulation time
CFL = 0.3;        % Desired CFL number 
dt = CFL * min(dx, dy); % Time step

% Advection speed function
a = @(x, y) [x -y]; % a(x) = (x1, -x2)

% Initial condition
u0 =@(x, y) zeros(size(x)); % Initial solution u(0, x) = 0

% Boundary conditions
bc_inflow_x1 = @(t, x2) 1 + x2.^2;    % Boundary condition at x1 = 1
bc_inflow_x2 = @(t, x1) 1 + 4*x1.^2;  % Boundary condition at x2 = 2

% Initialize grid and solution
x = linspace(1, 2, Nx);
y = linspace(1, 2, Ny);
[X, Y] = meshgrid(x, y);
u = u0(X, Y);

uR = 1 + (X .* Y).^2;

% Plot the function
surf(X, Y, uR);
% quiver(X, Y , X, -Y);

% Plot the final solution
figure;
surf(X, Y, u);
title('Final solution');
xlabel('x');
ylabel('y');
zlabel('u(x,y)');



% Time-stepping loop
t = 0;
while t < tFinal
    c =  a(X, Y);
    % Semi-Lagrangian update
    x_back = X - dt *c(:,1:Nx);  % Backward particle tracing in x-direction
    y_back = Y - dt * c(:,Nx+1:end);  % Backward particle tracing in y-direction

    % Apply boundary conditions at the inflow boundaries
    u(:, 1) = bc_inflow_x1(t, x);
    u(end,: ) = bc_inflow_x2(t, y);
    

    % Interpolation using LINEAR INTERPOLATION
    u_interp = interp2(X, Y, u, x_back, y_back, 'linear', 0);

    % Interpolation using RBF INTERPOLATION
%              pd=2;
%              ns=2*nchoosek(pd+2,2);
%              [ind,~] = knnsearch([X(:) Y(:)],[x_back(:) y_back(:)],'K',ns);
%              for j =1:numel(X)
%                  W = InterpRBF([X(ind(j,:))' Y(ind(j,:))'],[x_back(j) y_back(j)],1,pd);
%                  Z1(j) = sum(W.*u(ind(j,:))'); % approx. using RBFs
%              end
% %     
%              u_interp = reshape(Z1,[Nx Nx]);


    
    % Update solution
    u = u_interp;
    
%     surf(X, Y, u);
%     axis([1 Lx+1 1 Ly+1 -1 20]);
%     title(['Time = ', num2str(t)]);
%     xlabel('x');
%     ylabel('y');
%     zlabel('u');
%     drawnow;

    % Update time
    t = t + dt;


end


u(:, 1) = bc_inflow_x1(t, x);
u(end,: ) = bc_inflow_x2(t, y);
    
figure;
surf(X, Y, u);
title('Final solution');
xlabel('x');
ylabel('y');
zlabel('u(x,y)');
