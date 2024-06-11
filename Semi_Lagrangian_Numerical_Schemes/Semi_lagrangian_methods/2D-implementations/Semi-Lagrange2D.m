% Semi-Lagrangian Scheme for 2D Advection Equation
clear all;
clc;
points =[16, 32, 64, 90, 128];
count=1;
for i = 1:length(points)
    % Parameters
    Lx = 1;           % Length of the domain in x-direction
    Ly = 1;           % Length of the domain in y-direction
    Nx =  points(4);         % Number of grid points in x-direction
    Ny =  points(4);         % Number of grid points in y-direction
    dx = Lx / Nx;     % Grid spacing in x-direction
    dy = Ly / Ny;     % Grid spacing in y-direction
    cx = 1;            % Advection speed
    cy = 1; 

    % Time parameters
    tFinal = 5;     % Final simulation time
    CFL = 0.5;        % Desired CFL number 
%     dt = CFL * min(dx,dy)/sqrt(cx^2+cy^2); %RBF
    dt = CFL * dx;
%     dt =  CFL / (1/(dx) + 1/(dy)); % Quasicubic

    % Initial condition
    k = 10;
    x0 = 0.5;
    y0 = 0.5;
    u0 =@(x, y) exp(-k^2 * ((x-x0).^2 + (y-y0).^2));

    % Initialize grid and solution
    x = linspace(0, Lx, Nx);
    y = linspace(0, Ly, Ny);
    [X, Y] = meshgrid(x, y);
    u = u0(X, Y);
% 
    Nt = ceil(tFinal / dt);
    dt = tFinal/Nt;
%     Nt = floor(tFinal/dt);
%      dt = tFinal/Nt;
%     Vt = @(t)  pi/2 * sin(pi/2 * t + pi/2);
%      Vx = @(x,y)  0.2 + sin(x-1).^2  + cos(y-1).^2;



    %-----------------------


    % Time-stepping loop
    t = 0;
    while t < tFinal
        clear u1;

        % RK4 for time integration of the velocity field DEPENDING ON T
%         k1 = dt * Vt(t);
%         k2 = dt * Vt(t + 0.5 * dt);
%         k3 = dt * Vt(t + 0.5 * dt);
%         k4 = dt * Vt(t + dt);
% %        Semi-Lagrangian update
%          x_back = mod(X - (k1 + 2*k2 + 2*k3 + k4)/6, Lx);  % Backward particle tracing
%          y_back = mod(Y - (k1 + 2*k2 + 2*k3 + k4)/6, Ly);  % Backward particle tracing



        [X,Y,u] = periodic(X,Y,u); % impose periodic BCs
        [Xe,Ye,ue] = periodic_extension(X,Y,u); % periodic extension of data

        
        % Semi-Lagrangian update ---- (MCAMBIAR PERIODIC EXTENSION SEGUN METODO)
        x_back = mod(X - cx * dt, Lx);  % Backward particle tracing in x-direction
        y_back = mod(Y - cy * dt, Ly);  % Backward particle tracing in y-direction
        


%         %----------- RBF INTERPOLATION --------(MEJOR SIN PERIODIC EXTENSION)
%         pd=4;
%         ns=2*nchoosek(pd+2,2);
%         [ind,~] = knnsearch([Xe(:) Ye(:)],[x_back(:) y_back(:)],'K',ns);
%         for j =1:numel(Xe)
%             W = InterpRBF([Xe(ind(j,:))' Ye(ind(j,:))'],[x_back(j) y_back(j)],3,pd);
%             u1(j) = sum(W.*ue(ind(j,:))'); % approx. using RBFs
%         end
% 
% %           % u = reshape(u1,[Nx Nx])
%          u1 = reshape(u1,[Nx+2 Nx+2]);
%          u = u1;
%          u(:,end)=[];
%          u(:,1)=[];
%          u(end,:)=[];
%          u(1,:)=[];
%         %--------------------------------------

%           %----------- QUASICUBIC INTERPOLATION --------(PONER PERIODIC EXTENSION)
%         [ind,~] = knnsearch([Xe(:) Ye(:)],[x_back(:) y_back(:)],'K',4);
%       for j=1:numel(Xe)            
% %             u1(j) = quasicubic_interp_linear(x_back(j),y_back(j),Xe,Ye,ue,ind(j,:)); % approx using qusicubic interpolation
%             u1(j) = quasicubic_interp(x_back(j),y_back(j),Xe,Ye,ue,ind(j,:));
%       end
% 
% %       u = reshape(u1,[Nx Nx]);
%        u1 = reshape(u1,[Nx+2 Nx+2]);
%        u = u1;
%          u(:,end)=[];
%          u(:,1)=[];
%          u(end,:)=[];
%          u(1,:)=[];
% %          
         %--------------------------------------

        

        u = interp2(X, Y, u, x_back, y_back, 'spline', 0);  % 2D interpolation



        % Update time
        t = t + dt;
% 
%           % Plot the solution at each time step (optional)
%         surf(X, Y, u);
%         axis([0 Lx 0 Ly -1 1]);
%         title(['Time = ', num2str(t)]);
%         xlabel('x');
%         ylabel('y');
%         zlabel('u');
%         drawnow;

    end

%     ue = @(x, y) exp(-k^2 * ((x-x0).^2 + (y-y0).^2));
     Error = u-u0(X,Y);
    LINF(i) = norm(Error,inf);
    numErrorsL2(i)  = sqrt(sum(sum((Error).^2)) / (Nx^2));
    
         figure;
 
        surf(X, Y, u);
        axis([0 Lx 0 Ly -0.2 1.2]);
        xlabel('x');
        ylabel('y');
        zlabel('u');
        ax = gca;
        ax.FontSize = 15; 
        view([0 90]);
        shading interp;
        axis equal;
        colorbar;
         set(gca,'ytick',[0 0.2 0.4 0.6 0.8 1])
          set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1])
        drawnow;

   % L22(i) = sqrt(sum((ue(X,Y) - u).^2) / Nx*Ny);
    %sqrt(sum((u0 - u).^2) / Nx);

end

for i = 1:length(points)-1
    
    ordersConv(i) = log(numErrorsL2(i)/numErrorsL2(i+1))/log(points(i)/points(i+1));
end


% Plot the final solution

figure;
surf(X, Y, u);
title('Final solution');
xlabel('x');
ylabel('y');
zlabel('u(x,y)');


