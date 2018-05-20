function u = exact(x,y,t)
    %u = 500.0 + 500.0 * ((sin(pi*x))^2 + (sin(pi*y))^2)/2;
    u = 100*x^3*cos(2*pi*y)*cos(t)+20;
end