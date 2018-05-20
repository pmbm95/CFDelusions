function f = source(x,y,t)
    cd = -1;
    u = exact(x,y,t);
%     time = -1000*x^3*cos(2*pi*y)*sin(10*t+1);
%     conv_x = u*300*x^2*cos(2*pi*y)*cos(10*t+1);
%     conv_y = -u*200*pi*x^3*sin(2*pi*y)*cos(10*t+1);
%     diff_x = cd*600*x*cos(2*pi*y)*cos(10*t+1);
%     diff_y = -cd*400*pi^2*x^3*cos(2*pi*y)*cos(10*t+1);
%     f = time + conv_x + conv_y + diff_x + diff_y; 
 %time = 100*x^3*cos(2*pi*y)*exp(0.1*t)*0.1;
    conv_x = u*300*x^2*cos(2*pi*y);
    conv_y = -u*200*pi*x^3*sin(2*pi*y);
    diff_x = cd*600*x*cos(2*pi*y);
    diff_y = -cd*400*pi^2*x^3*cos(2*pi*y);
    f = conv_x + conv_y + diff_x + diff_y; 
end
