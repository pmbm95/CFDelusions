clc
clear all
close all

% Computation Time
% Vectorial
% Time Dependent
% Neumann Boundary Conditions
% Non Linear
% SOR/ Jacobi
% Residual

dominio_x = [-.5, .5];
dominio_y = [-.5, .5];
n_nodes_x = 60;
n_nodes_y = 60;

SOR = false;
w = 1.8;

time_dependent = false;
tmax = 15;
dt = 0.01;

dx = (dominio_x(2) - dominio_x(1) ) / n_nodes_x;
dy = (dominio_y(2) - dominio_y(1) ) / n_nodes_y;

centro_x = (dominio_x(1)+dx/2) : dx : (dominio_x(2)-dx/2);
centro_y = (dominio_y(1)+dy/2) : dy : (dominio_y(2)-dy/2);

[x_nodes, y_nodes] = meshgrid(centro_x, centro_y);

%plot(x_nodes, y_nodes, 'k.')

%--------- Coeff de difusao ----------
cd = -1;
if time_dependent
    Vx = arrayfun(@(u,v,t) exact(u,v,t), x_nodes, y_nodes, 0*ones(n_nodes_x,n_nodes_y));
    Vx = Vx';
    Vx = Vx(:);
    Vy = Vx;
else
    Vx = ones(n_nodes_x*n_nodes_y,1)*0;
    Vy = ones(n_nodes_x*n_nodes_y,1)*0;
end
S_boundary_n = false;
N_boundary_n = false;
W_boundary_n = false;
E_boundary_n = false;

S_n = 0;
N_n = 0;
W_n = 0;
E_n = 0;

if time_dependent
    U_old = arrayfun(@(x,y,t) exact(x,y,t), x_nodes, y_nodes, 0*ones(n_nodes_x, n_nodes_y));
    U_old = U_old';
    U_old = U_old(:);
    U_old2 = U_old;
else
    U_old = ones(n_nodes_x*n_nodes_y,1)*0;
end

for t=dt:dt:tmax
    
    A = sparse(n_nodes_y*n_nodes_x, n_nodes_y*n_nodes_x);
    B = sparse(n_nodes_y*n_nodes_x, 1);
    
    i = 1 ;% Vertical Counter
    j = 1; % Horizontal counter
    
%     if not(time_dependent)
%         t = 0;
%     end
    
    %Matrix Assembly
    for k=1:n_nodes_y*n_nodes_x
        % --- West Flux ---
        if j == 1 % Boundary
            %Diffusion
            A(k,k) = A(k,k) - 3.5*cd*(dy/dx);
            A(k,k+1) = A(k,k+1) + 0.5*cd*dy/dx;
            B(k) = B(k) - 3*cd*dy/dx*exact(centro_x(j)-(dx/2), centro_y(i),t);
            %Convection
            B(k) = B(k) + Vx(k)*dy*exact(centro_x(j)-(dx/2), centro_y(i),t);
        else % Center
            %Diffusion
            A(k,k) = A(k,k) - cd*(dy/dx);
            A(k,k-1) = A(k,k-1) + cd*(dy/dx);
            
            %Convection
            if Vx(k) > 0
                if j==2
                    A(k,k-1) = A(k,k-1) - 2*Vx(k)*dy;
                    B(k) = B(k) - 1*Vx(k)*dy*exact(centro_x(j-1)-(dx/2), centro_y(i),t);
                else
                    A(k,k-1)=A(k,k-1) - 3*(Vx(k)*dy*0.5);
                    A(k,k-2)=A(k,k-2) + 1*(Vx(k)*dy*0.5);
                end
            else %vento negativo
                if j== n_nodes_x
                    A(k,k) = A(k,k) - 2*Vx(k)*dy;
                    B(k) = B(k) - 1*Vx(k)*dy*exact(centro_x(j)+ (dx/2), centro_y(i),t);
                else
                    A(k,k)=A(k,k) - 3*(Vx(k)*dy*0.5);
                    A(k,k+1)=A(k,k+1) + 1*(Vx(k)*dy*0.5);
                end
            end
        end
        
        % --- East Flux ---
        if j == n_nodes_x % Boundary
            A(k,k)= A(k,k) - 3.5*cd*(dy/dx);
            A(k,k-1)= A(k,k-1) + 0.5*cd*(dy/dx);
            B(k)= B(k) - 3*cd*(dy/dx)*exact(centro_x(j)+(dx/2), centro_y(i),t);
            %Convection
            B(k) = B(k) - Vx(k)*dy*exact(centro_x(j)+(dx/2), centro_y(i),t);
        else % Center
            %Diffusion
            A(k,k)= A(k,k) - cd*(dy/dx);
            A(k,k+1)= A(k,k+1) + cd*(dy/dx);
            %Convection
            if Vx(k) > 0
                if j==1
                    A(k,k) = A(k,k) + 2*Vx(k)*dy;
                    B(k) = B(k) + 1*Vx(k)*dy*exact(centro_x(j)-(dx/2), centro_y(i),t);
                else
                    A(k,k)= A(k,k) + 3*(Vx(k)*dy*0.5);
                    A(k,k-1)= A(k,k-1) - 1*(Vx(k)*dy*0.5);
                end
            else %vento negativo
                if j==n_nodes_x-1
                    A(k,k+1)= A(k,k+1) + 2*(Vx(k)*dy);
                    B(k) = B(k) + 1*Vx(k)*dy*exact(centro_x(j+1)+(dx/2), centro_y(i),t);
                else
                    A(k,k+1)= A(k,k+1) + 3*(Vx(k)*dy*0.5);
                    A(k,k+2)= A(k,k+2) - 1*(Vx(k)*dy*0.5);
                end
            end
        end
        
        % --- South Flux ---
        if i == 1 % Boundary
            % Diffusion
            A(k,k)= A(k,k) - 3.5*cd*(dx/dy);
            A(k,k+n_nodes_x)= A(k,k+n_nodes_x) + 0.5*cd*(dx/dy);
            B(k)= B(k) - 3*cd*(dx/dy) * exact(centro_x(j), centro_y(i)-(dy/2),t);
            % Convection
            B(k) = B(k) + Vy(k)*dx*exact(centro_x(j), centro_y(i)-(dy/2),t);
        else % Center
            % Diffusion
            A(k,k) = A(k,k) - cd*(dx/dy);
            A(k,k-n_nodes_x) = A(k,k-n_nodes_x) + cd*(dx/dy);
            % Convection
            if Vy(k) > 0
                if i == 2
                    A(k,k-n_nodes_x) = A(k,k-n_nodes_x) - Vy(k)*dx*2;
                    B(k) = B(k) - Vy(k)*dx*exact(centro_x(j), centro_y(i-1)-(dy/2),t);
                else
                    A(k,k-n_nodes_x) = A(k,k-n_nodes_x) - 0.5*Vy(k)*dx*3;
                    A(k,k-2*n_nodes_x) = A(k,k-2*n_nodes_x) + 0.5*Vy(k)*dx;
                end
            else
                if i == n_nodes_y
                    A(k,k) = A(k,k) - Vy(k)*dx*2;
                    B(k) = B(k) - Vy(k)*dx*exact(centro_x(j), centro_y(i)+(dy/2),t);
                else
                    A(k,k) = A(k,k) - 0.5*Vy(k)*dx*3;
                    A(k,k+n_nodes_x) = A(k,k+n_nodes_x) + 0.5*Vy(k)*dx;
                end
            end
        end
        
        % --- North Flux ---
        if i == n_nodes_y % Boundary
            % Diffusion
            A(k,k) = A(k,k) - 3.5*cd*(dx/dy);
            A(k,k-n_nodes_x) = A(k,k-n_nodes_x) + 0.5*cd*dx/dy;
            B(k) = B(k) - 3*cd*(dx/dy)*exact(centro_x(j), centro_y(i)+(dy/2),t);
            % Convection
            B(k) = B(k) - Vy(k)*dx*exact(centro_x(j), centro_y(i)+(dy/2),t);
        else % Center
            % Diffusion
            A(k,k) = A(k,k) - cd*(dx/dy);
            A(k,k+n_nodes_x) = A(k,k+n_nodes_x) + cd*(dx/dy);
            % Convection
            if Vy(k) > 0
                if i == 1
                    A(k,k) = A(k,k) + Vy(k)*dx*2;
                    B(k) = B(k) + Vy(k)*dx*exact(centro_x(j), centro_y(i)-(dy/2),t);
                else
                    A(k,k) = A(k,k) + 0.5*Vy(k)*dx*3;
                    A(k,k-n_nodes_x) = A(k,k-n_nodes_x) - 0.5*Vy(k)*dx;
                end
            else
                if i == n_nodes_y-1
                    A(k,k+n_nodes_x) = A(k,k+n_nodes_x) + Vy(k)*dx*2;
                    B(k) = B(k) + Vy(k)*dx*exact(centro_x(j), centro_y(i+1)+(dy/2),t);
                else
                    A(k,k+n_nodes_x) = A(k,k+n_nodes_x) + 0.5*Vy(k)*dx*3;
                    A(k,k+2*n_nodes_x) = A(k,k+2*n_nodes_x) - 0.5*Vy(k)*dx;
                end
            end
        end
        
        
        % --- Source ---
        B(k) = B(k) + source(centro_x(j), centro_y(i),t)*dx*dy;
        
        
        % --- Time ---
        if time_dependent
            if t == dt
                A(k,k) = A(k,k) + dx*dy/dt;
                B(k) = B(k) + U_old(k)*dx*dy/dt;
            else
                A(k,k) = A(k,k) + 1.5*dx*dy/dt;
                B(k) = B(k) + 2*U_old(k)*dx*dy/dt - 0.5*U_old2(k)*dx*dy/dt;
            end
        end
          
        if j == n_nodes_x % Right Boundary
            j = 1;
            i = i+1;
        else
            j = j+1;
        end
    end
    
    if SOR
        n_iter_limit = 300;
        n_iter = 1;
        residual = 1;
        Low = -tril(A,-1);
        Up = -triu(A,1);
        D = diag(A);
        u_old = U_old;
        u = zeros(n_nodes_x*n_nodes_y,1);
        while residual >= 1e-3
            for i=1:n_nodes_x*n_nodes_y
                u(i) = D(i)^-1 * ( Low(i,:)*u + Up(i,:)*u_old + B(i));
                u(i) = u_old(i) + w*(u(i)-u_old(i));
            end
            %residual = sum((u-u_old).^2)^0.5;
            residual = max(abs(u-u_old));
            u_old = u;
            
            display = ['Iteration: ', num2str(n_iter)];
            disp(display)
            display = ['Residual: ', num2str(residual)];
            disp(display)
            
            if n_iter >= n_iter_limit
                break
            end
            n_iter = n_iter + 1;
        end
        U = u;
    else
        U = A\B;
    end
    
    U_old = U;
    Vx=U;
    Vy=U;
    if time_dependent
        U_old2 = U_old;
    end
    
    U = vec2mat(U,n_nodes_x);
    
    U_exact = arrayfun(@(u,v,t) exact(u,v,t), x_nodes, y_nodes, t*ones(n_nodes_x,n_nodes_y));
    
    Error = abs((U_exact - U));
    
    figure(1)
    h = pcolor(x_nodes, y_nodes, U);
    title(['U at time: ', num2str(t)]);
    %set(h, 'EdgeColor', 'none');
    colormap(jet);
    colorbar;
    %caxis([5, 35]);
    
    figure(2)
    h = pcolor(x_nodes, y_nodes, Error);
    %set(h, 'EdgeColor', 'none');
    title(['Error at time: ', num2str(t)]);
    colormap(jet);
    colorbar;
    %caxis([0, .12]);
    
    %max(max(Error))
    
    pause(0.05);
    
   
end