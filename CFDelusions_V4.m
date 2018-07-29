clc
clear all
close all

% Computation Time

% Neumann Boundary Conditions
% Residual

dominio_x = [-.5, .5];
dominio_y = [-.5, .5];

n_nodes_x = 100;
n_nodes_y = 100;

SOR = true;
w = 1.5;

time_dependent = true;
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

Vx = x_nodes*0;
Vy = x_nodes*0;
if time_dependent
    for i=1:n_nodes_y
        for j=1:n_nodes_x   
        a = exact(x_nodes(i,j), y_nodes(i,j), 0);
        Vx(i,j) = a(1);
        Vy(i,j) = a(2);
        end
    end
    Vx = Vx';
    Vx = Vx(:);
    Vy = Vy';
    Vy = Vy(:);
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
    U_old = Vx;
    V_old = Vy;
    U_old = U_old';
    U_old = U_old(:);
    U_old2 = U_old;
    V_old = V_old';
    V_old = V_old(:);
    V_old2 = V_old;
else
    U_old = ones(n_nodes_x*n_nodes_y,1)*0;
    V_old = ones(n_nodes_x*n_nodes_y,1)*0;
end

for t=dt:dt:tmax
    
    
    A = sparse(n_nodes_y*n_nodes_x, n_nodes_y*n_nodes_x);
    Bx = sparse(n_nodes_y*n_nodes_x, 1);
    Ay = sparse(n_nodes_y*n_nodes_x, n_nodes_y*n_nodes_x);
    By = sparse(n_nodes_y*n_nodes_x, 1);
    i = 1 ;% Vertical Counter
    j = 1; % Horizontal counter
    
    if not(time_dependent)
        t = tmax;
    end
    
    %Matrix Assembly PRECOMPUTATION
    for k=1:n_nodes_y*n_nodes_x
        % --- West Flux ---
        if j == 1 % Boundary
            %Diffusion
            A(k,k) = A(k,k) - 3.5*cd*(dy/dx);
            A(k,k+1) = A(k,k+1) + 0.5*cd*dy/dx;
            Bx(k) = Bx(k) - 3*cd*dy/dx*exact(centro_x(j)-(dx/2), centro_y(i),t)*[1;0];
            By(k) = By(k) - 3*cd*dy/dx*exact(centro_x(j)-(dx/2), centro_y(i),t)*[0;1];
        else % Center
            %Diffusion
            A(k,k) = A(k,k) - cd*(dy/dx);
            A(k,k-1) = A(k,k-1) + cd*(dy/dx);
        end
        
        % --- East Flux ---
        if j == n_nodes_x % Boundary
            A(k,k)= A(k,k) - 3.5*cd*(dy/dx);
            A(k,k-1)= A(k,k-1) + 0.5*cd*(dy/dx);
            Bx(k)= Bx(k) - 3*cd*(dy/dx)*exact(centro_x(j)+(dx/2), centro_y(i),t)*[1;0];
            By(k)= By(k) - 3*cd*(dy/dx)*exact(centro_x(j)+(dx/2), centro_y(i),t)*[0;1];
        else % Center
            %Diffusion
            A(k,k)= A(k,k) - cd*(dy/dx);
            A(k,k+1)= A(k,k+1) + cd*(dy/dx);
        end
        
        % --- South Flux ---
        if i == 1 % Boundary
            % Diffusion
            A(k,k)= A(k,k) - 3.5*cd*(dx/dy);
            A(k,k+n_nodes_x)= A(k,k+n_nodes_x) + 0.5*cd*(dx/dy);
            Bx(k)= Bx(k) - 3*cd*(dx/dy) * exact(centro_x(j), centro_y(i)-(dy/2),t)*[1;0];
            By(k)= By(k) - 3*cd*(dx/dy) * exact(centro_x(j), centro_y(i)-(dy/2),t)*[0;1];
        else % Center
            % Diffusion
            A(k,k) = A(k,k) - cd*(dx/dy);
            A(k,k-n_nodes_x) = A(k,k-n_nodes_x) + cd*(dx/dy);
        end
        
        % --- North Flux ---
        if i == n_nodes_y % Boundary
            % Diffusion
            A(k,k) = A(k,k) - 3.5*cd*(dx/dy);
            A(k,k-n_nodes_x) = A(k,k-n_nodes_x) + 0.5*cd*dx/dy;
            Bx(k) = Bx(k) - 3*cd*(dx/dy)*exact(centro_x(j), centro_y(i)+(dy/2),t)*[1;0];
            By(k) = By(k) - 3*cd*(dx/dy)*exact(centro_x(j), centro_y(i)+(dy/2),t)*[0;1];
        else % Center
            % Diffusion
            A(k,k) = A(k,k) - cd*(dx/dy);
            A(k,k+n_nodes_x) = A(k,k+n_nodes_x) + cd*(dx/dy);
        end
        
        % --- Source ---
        Bx(k) = Bx(k) + source(centro_x(j), centro_y(i),t)*[1;0]*dx*dy;
        By(k) = By(k) + source(centro_x(j), centro_y(i),t)*[0;1]*dx*dy;
        
        % --- Time ---
        if time_dependent
            if t == dt
                A(k,k) = A(k,k) + dx*dy/dt;
                Bx(k) = Bx(k) + U_old(k)*dx*dy/dt;
                By(k) = By(k) + V_old(k)*dx*dy/dt;
            else
                A(k,k) = A(k,k) + 1.5*dx*dy/dt;
                Bx(k) = Bx(k) + 2*U_old(k)*dx*dy/dt - 0.5*U_old2(k)*dx*dy/dt;
                By(k) = By(k) + 2*V_old(k)*dx*dy/dt - 0.5*V_old2(k)*dx*dy/dt;
            end
        end
        
        if j == n_nodes_x % Right Boundary
            j = 1;
            i = i+1;
        else
            j = j+1;
        end
    end
    
    A_pre = A;
    Bx_pre = Bx;
    By_pre = By;
    n_iter_limit = 300;
    n_iter = 1;
    residual = 1;
    u_old = U_old;
    v_old = V_old;
    while residual >= 1e-3
        
        A = A_pre;
        Bx = Bx_pre;
        By = By_pre;
        %Matrix Assembly NON LINEAR
        i = 1 ;% Vertical Counter
        j = 1; % Horizontal counter
        for k=1:n_nodes_y*n_nodes_x
            % --- West Flux ---
            if j == 1 % Boundary
                %Convection
                Bx(k) = Bx(k) + Vx(k)*dy*exact(centro_x(j)-(dx/2), centro_y(i),t)*[1;0];
                By(k) = By(k) + Vx(k)*dy*exact(centro_x(j)-(dx/2), centro_y(i),t)*[0;1];
            else % Center
                %Convection
                if Vx(k) > 0
                    if j==2
                        A(k,k-1) = A(k,k-1) - 2*Vx(k)*dy;
                        Bx(k) = Bx(k) - 1*Vx(k)*dy*exact(centro_x(j-1)-(dx/2), centro_y(i),t)*[1;0];
                        By(k) = By(k) - 1*Vx(k)*dy*exact(centro_x(j-1)-(dx/2), centro_y(i),t)*[0;1];
                    else
                        A(k,k-1)=A(k,k-1) - 3*(Vx(k)*dy*0.5);
                        A(k,k-2)=A(k,k-2) + 1*(Vx(k)*dy*0.5);
                    end
                else %vento negativo
                    if j== n_nodes_x
                        A(k,k) = A(k,k) - 2*Vx(k)*dy;
                        Bx(k) = Bx(k) - 1*Vx(k)*dy*exact(centro_x(j)+ (dx/2), centro_y(i),t)*[1;0];
                        By(k) = By(k) - 1*Vx(k)*dy*exact(centro_x(j)+ (dx/2), centro_y(i),t)*[0;1];
                    else
                        A(k,k)=A(k,k) - 3*(Vx(k)*dy*0.5);
                        A(k,k+1)=A(k,k+1) + 1*(Vx(k)*dy*0.5);
                    end
                end
            end
            
            % --- East Flux ---
            if j == n_nodes_x % Boundary
                %Convection
                Bx(k) = Bx(k) - Vx(k)*dy*exact(centro_x(j)+(dx/2), centro_y(i),t)*[1;0];
                By(k) = By(k) - Vx(k)*dy*exact(centro_x(j)+(dx/2), centro_y(i),t)*[0;1];
            else % Center
                %Convection
                if Vx(k) > 0
                    if j==1
                        A(k,k) = A(k,k) + 2*Vx(k)*dy;
                        Bx(k) = Bx(k) + 1*Vx(k)*dy*exact(centro_x(j)-(dx/2), centro_y(i),t)*[1;0];
                        By(k) = By(k) + 1*Vx(k)*dy*exact(centro_x(j)-(dx/2), centro_y(i),t)*[0;1];
                    else
                        A(k,k)= A(k,k) + 3*(Vx(k)*dy*0.5);
                        A(k,k-1)= A(k,k-1) - 1*(Vx(k)*dy*0.5);
                    end
                else %vento negativo
                    if j==n_nodes_x-1
                        A(k,k+1)= A(k,k+1) + 2*(Vx(k)*dy);
                        Bx(k) = Bx(k) + 1*Vx(k)*dy*exact(centro_x(j+1)+(dx/2), centro_y(i),t)*[1;0];
                        By(k) = By(k) + 1*Vx(k)*dy*exact(centro_x(j+1)+(dx/2), centro_y(i),t)*[0;1];
                    else
                        A(k,k+1)= A(k,k+1) + 3*(Vx(k)*dy*0.5);
                        A(k,k+2)= A(k,k+2) - 1*(Vx(k)*dy*0.5);
                    end
                end
            end
            
            % --- South Flux ---
            if i == 1 % Boundary
                % Convection
                Bx(k) = Bx(k) + Vy(k)*dx*exact(centro_x(j), centro_y(i)-(dy/2),t)*[1;0];
                By(k) = By(k) + Vy(k)*dx*exact(centro_x(j), centro_y(i)-(dy/2),t)*[0;1];
            else % Center
                % Convection
                if Vy(k) > 0
                    if i == 2
                        A(k,k-n_nodes_x) = A(k,k-n_nodes_x) - Vy(k)*dx*2;
                        Bx(k) = Bx(k) - Vy(k)*dx*exact(centro_x(j), centro_y(i-1)-(dy/2),t)*[1;0];
                        By(k) = By(k) - Vy(k)*dx*exact(centro_x(j), centro_y(i-1)-(dy/2),t)*[0;1];
                    else
                        A(k,k-n_nodes_x) = A(k,k-n_nodes_x) - 0.5*Vy(k)*dx*3;
                        A(k,k-2*n_nodes_x) = A(k,k-2*n_nodes_x) + 0.5*Vy(k)*dx;
                    end
                else
                    if i == n_nodes_y
                        A(k,k) = A(k,k) - Vy(k)*dx*2;
                        Bx(k) = Bx(k) - Vy(k)*dx*exact(centro_x(j), centro_y(i)+(dy/2),t)*[1;0];
                        By(k) = By(k) - Vy(k)*dx*exact(centro_x(j), centro_y(i)+(dy/2),t)*[0;1];
                    else
                        A(k,k) = A(k,k) - 0.5*Vy(k)*dx*3;
                        A(k,k+n_nodes_x) = A(k,k+n_nodes_x) + 0.5*Vy(k)*dx;
                    end
                end
            end
            
            % --- North Flux ---
            if i == n_nodes_y % Boundary
                % Convection
                Bx(k) = Bx(k) - Vy(k)*dx*exact(centro_x(j), centro_y(i)+(dy/2),t)*[1;0];
                By(k) = By(k) - Vy(k)*dx*exact(centro_x(j), centro_y(i)+(dy/2),t)*[0;1];
            else % Center
                % Convection
                if Vy(k) > 0
                    if i == 1
                        A(k,k) = A(k,k) + Vy(k)*dx*2;
                        Bx(k) = Bx(k) + Vy(k)*dx*exact(centro_x(j), centro_y(i)-(dy/2),t)*[1;0];
                        By(k) = By(k) + Vy(k)*dx*exact(centro_x(j), centro_y(i)-(dy/2),t)*[0;1];
                    else
                        A(k,k) = A(k,k) + 0.5*Vy(k)*dx*3;
                        A(k,k-n_nodes_x) = A(k,k-n_nodes_x) - 0.5*Vy(k)*dx;
                    end
                else
                    if i == n_nodes_y-1
                        A(k,k+n_nodes_x) = A(k,k+n_nodes_x) + Vy(k)*dx*2;
                        Bx(k) = Bx(k) + Vy(k)*dx*exact(centro_x(j), centro_y(i+1)+(dy/2),t)*[1;0];
                        By(k) = By(k) + Vy(k)*dx*exact(centro_x(j), centro_y(i+1)+(dy/2),t)*[0;1];
                    else
                        A(k,k+n_nodes_x) = A(k,k+n_nodes_x) + 0.5*Vy(k)*dx*3;
                        A(k,k+2*n_nodes_x) = A(k,k+2*n_nodes_x) - 0.5*Vy(k)*dx;
                    end
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
            Low = -tril(A,-1);
            Up = -triu(A,1);
            D = diag(A);
            u = zeros(n_nodes_x*n_nodes_y,1);
            v = zeros(n_nodes_x*n_nodes_y,1);
            for i=1:n_nodes_x*n_nodes_y
                u(i) = D(i)^-1 * ( Low(i,:)*u + Up(i,:)*u_old + Bx(i));
                u(i) = u_old(i) + w*(u(i)-u_old(i));
                v(i) = D(i)^-1 * ( Low(i,:)*v + Up(i,:)*v_old + By(i));
                v(i) = v_old(i) + w*(v(i)-v_old(i));
            end
            U = u;
            V = v;
        else
            U = A\Bx;
            V = A\By;
        end
        
        %residual = sum((U-u_old).^2)^0.5;
        residualx = max(abs(U-u_old));
        residualy = max(abs(V-v_old));
        residual = max([residualx residualy]);
        u_old = U;
        v_old = V;
        Vx=U;
        Vy=V;
        
        display = ['Iteration: ', num2str(n_iter)];
        disp(display)
        display = ['Residual: ', num2str(residual)];
        disp(display)
        
        if n_iter >= n_iter_limit
            break
        end
        n_iter = n_iter + 1;
    end
   
    
    if time_dependent
        U_old2 = U_old;
        V_old2 = V_old;
    end
    U_old = U;
    V_old = V;
    
    
    U = vec2mat(U,n_nodes_x);
    V = vec2mat(V,n_nodes_x);
    
    U_exact = U*0;
    V_exact = V*0;
    for i=1:n_nodes_y
        for j=1:n_nodes_x   
        a = exact(x_nodes(i,j), y_nodes(i,j), t);
        U_exact(i,j) = a(1);
        V_exact(i,j) = a(2);
        end
    end
    Errorx = abs((U_exact - U));
    Errory = abs((V_exact - V));
    
    figure(1)
    h = pcolor(x_nodes, y_nodes, U);
    title(['U at time: ', num2str(t)]);
    %set(h, 'EdgeColor', 'none');
    colormap(jet);
    colorbar;
    %caxis([5, 35]);
    
    figure(2)
    h = pcolor(x_nodes, y_nodes, Errorx);
    %set(h, 'EdgeColor', 'none');
    title(['Error X at time: ', num2str(t)]);
    colormap(jet);
    colorbar;
    %caxis([0, .12]);
    
    figure(3)
    h = pcolor(x_nodes, y_nodes, V);
    title(['V at time: ', num2str(t)]);
    %set(h, 'EdgeColor', 'none');
    colormap(jet);
    colorbar;
    %caxis([5, 35]);
    
    figure(4)
    h = pcolor(x_nodes, y_nodes, Errory);
    %set(h, 'EdgeColor', 'none');
    title(['Error Y at time: ', num2str(t)]);
    colormap(jet);
    colorbar;
    %caxis([0, .12]);
    
    %max(max(Error))
    
    pause(0.0001);
    
    
end