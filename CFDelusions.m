dominio_x = [-0.5 0.5];
dominio_y = [-0.5 0.5];
n_nodes_x = 20;
n_nodes_y = 10;

dx = (dominio_x(2) - dominio_x(1) ) / n_nodes_x;
dy = (dominio_y(2) - dominio_y(1) ) / n_nodes_y;

centro_x = (dominio_x(1)+dx/2) : dx : (dominio_x(2)-dx/2);
centro_y = (dominio_y(1)+dy/2) : dy : (dominio_y(2)-dy/2);

[x_nodes y_nodes] = meshgrid(centro_x, centro_y);

plot(x_nodes, y_nodes, 'k.')

A = zeros(n_nodes_y*n_nodes_x);
B = zeros(n_nodes_y*n_nodes_x, 1);
%---------Ceoff de difusao para porjetos futuros----------
cd = 1;

for k=1:n_nodes_y*n_nodes_x
    
    A(k,k) = A(k,k) - 2 *cd*(dx/dy + dy/dx);
    A(k,k-1) = A(k,k-1) + cd*dy/dx;
    A(k,k+1) = A(k,k+1) + cd*dy/dx;
    A(k,k+n_nodes_x) = A(k,k+n_nodes_x) + cd*dx/dy;
    A(k,k-n_nodes_x) = A(k,k-n_nodes_x) + cd*dx/dy;
    
    
end

