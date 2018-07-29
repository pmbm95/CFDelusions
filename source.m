function f = source(x,y,t)
    cd = -1;
    a = exact(x,y,t);
    vx = a(1);
    vy = a(2);
%     time = -1000*x^3*cos(2*pi*y)*sin(10*t+1);
%     conv_x = u*300*x^2*cos(2*pi*y)*cos(10*t+1);
%     conv_y = -u*200*pi*x^3*sin(2*pi*y)*cos(10*t+1);
%     diff_x = cd*600*x*cos(2*pi*y)*cos(10*t+1);
%     diff_y = -cd*400*pi^2*x^3*cos(2*pi*y)*cos(10*t+1);
%     f = time + conv_x + conv_y + diff_x + diff_y; 
 %time = 100*x^3*cos(2*pi*y)*exp(0.1*t)*0.1;
    conv_x = vx*300*x^2*cos(2*pi*y);
    conv_y = -vy*200*pi*x^3*sin(2*pi*y);
    diff_x = cd*600*x*cos(2*pi*y);
    diff_y = -cd*400*pi^2*x^3*cos(2*pi*y);
    fx = conv_x + conv_y + diff_x + diff_y + 10; 

    conv_y = vy*300*y^2*cos(2*pi*x);
    conv_x = -vx*200*pi*y^3*sin(2*pi*x);
    diff_y = cd*600*y*cos(2*pi*x);
    diff_x = -cd*400*pi^2*y^3*cos(2*pi*x);
    fy = conv_x + conv_y + diff_x + diff_y +20; 
    f = [fx fy];
end
