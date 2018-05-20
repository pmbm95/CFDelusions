function u = exact(x,y,t)
    %u = 500.0 + 500.0 * ((sin(pi*x))^2 + (sin(pi*y))^2)/2;
    u = 100*x^3*cos(2*pi*y)*exp(-0.1*t)+20;
    %u = 50* x^2 *  y^2 + 5;
    %u = 500*x.^3 + 500*y.^3 + 500*x.^2 + 500*sin((y)*pi);
    %u = x^2 + y^2 + 20;
end