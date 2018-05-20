function f = source(x,y,t)
    cd = -1;
    u = exact(x,y,t);
    %f = ((300*x^2*Vx-100*x^3*4*pi^2*cd+600*x*cd)*cos(2*pi*y)-100*x^3*sin(2*pi*y)*2*pi*Vy)*exp(-0.1*t)-10*x^3*cos(2*pi*y)*exp(-0.1*t);
    time = -100*x^3*cos(2*pi*y)*sin(t);
    conv_x = u*300*x^2*cos(2*pi*y)*cos(t);
    conv_y = -u*200*pi*x^3*sin(2*pi*y)*cos(t);
    diff_x = cd*600*x*cos(2*pi*y)*cos(t);
    diff_y = -cd*400*pi^2*x^3*cos(2*pi*y)*cos(t);
    f = time + conv_x + conv_y + diff_x + diff_y; 
end
