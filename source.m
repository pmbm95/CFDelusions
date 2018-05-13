function f = source(x,y)
    Vx = 100;
    Vy = 100;
    cd = -1;
    %f =  500.0 * pi^2 * (cos(2*pi*x) + cos(2*pi*y));
    f = (300*x^2*Vx-100*x^3*4*pi^2*cd+600*x*cd)*cos(2*pi*y)-100*x^3*sin(2*pi*y)*2*pi*Vy;
    %f = 100*x*y*(y*Vx+x*Vy) + 100*cd*(x^2+y^2);
    %f = 500*6*x + 500*6*y + 500*2 - 500*pi*pi*sin((y)*pi);
    %f = 4+ 0*x + 0*y;
end
