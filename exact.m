function u = exact(x,y,t)
    %u = 500.0 + 500.0 * ((sin(pi*x))^2 + (sin(pi*y))^2)/2;
    vx = 100*x^3*cos(2*pi*y)+10;
    vy = 100*y^3*cos(2*pi*x)+20;
    u = vx;
end
