function u = exact(x,y)
    u = 500.0 + 500.0 * ((sin(pi*x))^2 + (sin(pi*y))^2)/2;
    %u = 500*x.^3 + 500*y.^3 + 500*x.^2; %+ 10*sin((y)*pi);
end