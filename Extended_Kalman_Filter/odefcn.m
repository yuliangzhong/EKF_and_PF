% In this file, we define the ODE function
% y_dot = q(y,v=0,t)
% y = [px, py, sx, sy, phi, rou, b]';
function dydt = odefcn(t,y,u,const)
Cdh = const.dragCoefficientHydr; % C_d,h
Cda = const.dragCoefficientAir;  % C_d,a
Cr = const.rudderCoefficient; % C_r
Cw = const.windVel; % C_w
px = y(1);
py = y(2);
sx = y(3);
sy = y(4);
phi = y(5);
rou = y(6);
b = y(7);

dydt = zeros(1,7);
dydt(1) = sx;
dydt(2) = sy;
dydt(3) = cos(phi)*(tanh(u(1))-Cdh*(sx^2+sy^2))...
          -Cda*(sx-Cw*cos(rou))*sqrt((sx-Cw*cos(rou))^2 + (sy-Cw*sin(rou))^2);
dydt(4) = sin(phi)*(tanh(u(1))-Cdh*(sx^2+sy^2))...
          -Cda*(sx-Cw*sin(rou))*sqrt((sx-Cw*cos(rou))^2 + (sy-Cw*sin(rou))^2);
dydt(5) = Cr*u(2);
dydt(6) = 0;
dydt(7) = 0;