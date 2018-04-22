clc
clear all

dominio_x = [-0.5 0.5];
dominio_y = [-0.5 0.5];
n_nodes_x = 20;
n_nodes_y = 10;

dx = (dominio_x(2) - dominio_x(1) ) / n_nodes_x;
dy = (dominio_y(2) - dominio_y(1) ) / n_nodes_y;

centro_x = (dominio_x(1)+dx/2) : dx : (dominio_x(2)-dx/2);
centro_y = (dominio_y(1)+dy/2) : dy : (dominio_y(2)-dy/2);

[x_nodes, y_nodes] = meshgrid(centro_x, centro_y);

%plot(x_nodes, y_nodes, 'k.')

A = zeros(n_nodes_y*n_nodes_x);
B = zeros(n_nodes_y*n_nodes_x, 1);
%---------Ceoff de difusao para projetos futuros----------
cd = 1;

i = 1 ;% Vertical Counter
j = 1; % Horizontal counter


for k=1:n_nodes_y*n_nodes_x
    
    % --- West Flux ---
    if j == 1 % Boundary
        A(k,k) = A(k,k) - 3.5*cd*dy/dx;
        A(k,k+1) = A(k,k+1) + 0.5*cd*dy/dx;
        B(k) = B(k)- 3*cd*dy/dx*exact(centro_x(j)-dx/2, centro_y(i));
    else % Center
        A(k,k) = A(k,k) - cd*dy/dx;
        A(k,k-1) = A(k,k-1) + cd*dy/dx;        
    end
    
    if j == n_nodes_x % East Flux  
       %e boundary
       A(k,k)= A(k,k) - 3.5*cd*(dy/dx);
       A(k,k-1)= A(k,k-1) + 0.5*cd*(dy/dx);
       B(k)= B(k) - 3*cd*(dy/dx)*exact(centro_x(j)+(dx/2), centro_y(i));
    else
        % e normal
        A(k,k)= A(k,k) - cd*(dy/dx);
        A(k,k+1)= A(k,k+1) + cd*(dy/dx);
    end
    
    if i == 1 % South Flux
      %e boundary
      
    else
        % e normal
        A(k,k) = A(k,k) - cd*(dx/dy);
        A(k,k-n_nodes_x) = A(k,k-n_nodes_x) + cd*(dx/dy);
    end
    
    if i == n_nodes_y % North Flux
    %e boundary
    else
        % e normal
    end
    
    % Middle
    A(k,k) = A(k,k) - 2 *cd*(dx/dy + dy/dx);
    A(k,k-1) = A(k,k-1) + cd*dy/dx;
    A(k,k+1) = A(k,k+1) + cd*dy/dx;
    A(k,k+n_nodes_x) = A(k,k+n_nodes_x) + cd*dx/dy;
    A(k,k-n_nodes_x) = A(k,k-n_nodes_x) + cd*dx/dy;
    
    
    B(k) = B(k) - source(centro_x(j), centro_y(i));
    
    if j == n_nodes_x % Right Boundary
        j = 1;
        i = i+1;
    else
        j = j+1;
    end
end

fprintf('lol')
