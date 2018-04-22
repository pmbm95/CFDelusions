clc
clear all
close all

dominio_x = [-0.5 0.5];
dominio_y = [-0.5 0.5];
n_nodes_x = 100;
n_nodes_y = 100;

dx = (dominio_x(2) - dominio_x(1) ) / n_nodes_x;
dy = (dominio_y(2) - dominio_y(1) ) / n_nodes_y;

centro_x = (dominio_x(1)+dx/2) : dx : (dominio_x(2)-dx/2);
centro_y = (dominio_y(1)+dy/2) : dy : (dominio_y(2)-dy/2);

[x_nodes, y_nodes] = meshgrid(centro_x, centro_y);

%plot(x_nodes, y_nodes, 'k.')

A = sparse(n_nodes_y*n_nodes_x, n_nodes_y*n_nodes_x);
B = sparse(n_nodes_y*n_nodes_x, 1);

%--------- Coeff de difusao ----------
cd = 1;

i = 1 ;% Vertical Counter
j = 1; % Horizontal counter


S_boundary_n = false;
N_boundary_n = false;
W_boundary_n = false;
E_boundary_n = false;

S_n = 0;
N_n = 0;
W_n = 0;
E_n = 0;

for k=1:n_nodes_y*n_nodes_x
    
    % --- West Flux ---
    if j == 1 % Boundary
        A(k,k) = A(k,k) - 3.5*cd*dy/dx;
        A(k,k+1) = A(k,k+1) + 0.5*cd*dy/dx;
        B(k) = B(k)- 3*cd*dy/dx*exact(centro_x(j)-(dx/2), centro_y(i));
        %centro_x(j)-(dx/2)
        %centro_y(i)
    else % Center
        A(k,k) = A(k,k) - cd*dy/dx;
        A(k,k-1) = A(k,k-1) + cd*dy/dx;        
    end
    
    % --- East Flux --- 
    if j == n_nodes_x % Boundary
       A(k,k)= A(k,k) - 3.5*cd*(dy/dx);
       A(k,k-1)= A(k,k-1) + 0.5*cd*(dy/dx);
       B(k)= B(k) - 3*cd*(dy/dx)*exact(centro_x(j)+(dx/2), centro_y(i));
    else % Center
        A(k,k)= A(k,k) - cd*(dy/dx);
        A(k,k+1)= A(k,k+1) + cd*(dy/dx);
    end
    
    % --- South Flux ---
    if i == 1 % Boundary
        A(k,k)= A(k,k) - 3.5*cd*(dy/dx);
        A(k,k+n_nodes_x)= A(k,k+n_nodes_x) + 0.5*cd*(dy/dx);
        B(k)= B(k) - 3*cd*(dy/dx) * exact(centro_x(j), centro_y(i)-(dy/2));
    else % Center
        A(k,k) = A(k,k) - cd*(dx/dy);
        A(k,k-n_nodes_x) = A(k,k-n_nodes_x) + cd*(dx/dy);
    end
    
    % --- North Flux ---
    if i == n_nodes_y % Boundary
        A(k,k) = A(k,k) - 3.5*cd*dx/dy;
        A(k,k-n_nodes_x) = A(k,k-n_nodes_x) + 0.5*cd*dy/dy;
        B(k) = B(k) - 3*cd*dx/dy*exact(centro_x(j), centro_y(i)+(dy/2));
    else % Center
        A(k,k) = A(k,k) - cd*dx/dy;
        A(k,k+n_nodes_x) = A(k,k+n_nodes_x) + cd*dx/dy;
    end
        % --- Source---

    B(k) = B(k) - source(centro_x(j), centro_y(i));
    
    if j == n_nodes_x % Right Boundary
        j = 1;
        i = i+1;
    else
        j = j+1;
    end
end



U = A\B;

U = vec2mat(U,n_nodes_x);
pcolor(x_nodes, y_nodes, U)
colorbar;
U_exact = exact(x_nodes, y_nodes);

figure
Error = U_exact - U;
pcolor(x_nodes, y_nodes, U_exact);
colorbar;


%fprintf('lol')
