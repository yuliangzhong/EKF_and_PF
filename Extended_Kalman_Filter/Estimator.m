function [posEst,linVelEst,oriEst,windEst,driftEst,...
          posVar,linVelVar,oriVar,windVar,driftVar,estState] = ...
    Estimator(estState,actuate,sense,tm,estConst)
% [posEst,linVelEst,oriEst,windEst,driftEst,...
%  posVar,linVelVar,oriVar,windVar,driftVar,estState] = 
% Estimator(estState,actuate,sense,tm,estConst)
%
% The estimator.
%
% The function is initialized for tm == 0, otherwise the estimator does an
% iteration step (compute estimates for the time step k).
%
% Inputs:
%   estState        previous estimator state (time step k-1)
%                   May be defined by the user (for example as a struct).
%   actuate         control input u(k-1), [1x2]-vector
%                   actuate(1): u_t, thrust command
%                   actuate(2): u_r, rudder command
%   sense           sensor measurements z(k), [1x5]-vector, INF entry if no
%                   measurement
%                   sense(1): z_a, distance measurement a
%                   sense(2): z_b, distance measurement b
%                   sense(3): z_c, distance measurement c
%                   sense(4): z_g, gyro measurement
%                   sense(5): z_n, compass measurement
%   tm              time t_k, scalar
%                   If tm==0 initialization, otherwise estimator
%                   iteration step.
%   estConst        estimator constants (as in EstimatorConst.m)
%
% Outputs:
%   posEst          position estimate (time step k), [1x2]-vector
%                   posEst(1): p_x position estimate
%                   posEst(2): p_y position estimate
%   linVelEst       velocity estimate (time step k), [1x2]-vector
%                   linVelEst(1): s_x velocity estimate
%                   linVelEst(2): s_y velocity estimate
%   oriEst          orientation estimate (time step k), scalar
%   windEst         wind direction estimate (time step k), scalar
%   driftEst        estimate of the gyro drift b (time step k), scalar
%   posVar          variance of position estimate (time step k), [1x2]-vector
%                   posVar(1): x position variance
%                   posVar(2): y position variance
%   linVelVar       variance of velocity estimate (time step k), [1x2]-vector
%                   linVelVar(1): x velocity variance
%                   linVelVar(2): y velocity variance
%   oriVar          variance of orientation estimate (time step k), scalar
%   windVar         variance of wind direction estimate(time step k), scalar
%   driftVar        variance of gyro drift estimate (time step k), scalar
%   estState        current estimator state (time step k)
%                   Will be input to this function at the next call.
%
%
% Class:
% Recursive Estimation
% Spring 2021
% Programming Exercise 1
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Matthias Hofer, Carlo Sferrazza
% hofermat@ethz.ch
% csferrazza@ethz.ch

%% Initialization
if (tm == 0)
    % Do the initialization of your estimator here!
    R = estConst.StartRadiusBound;
    p_bar = estConst.RotationStartBound;
    r_bar = estConst.WindAngleStartBound;
    
    % initial state mean
    posEst = [0,0]; % 1x2 matrix
    linVelEst = [0,0]; % 1x2 matrix
    oriEst = 0; % 1x1 matrix
    windEst = 0; % 1x1 matrix
    driftEst = 0; % 1x1 matrix
    
    % initial state variance
    posVar = [pi*(R^4)/4,pi*(R^4)/4]; % 1x2 matrix
    linVelVar = [0,0]; % 1x2 matrix
    oriVar = (p_bar^2)/3; % 1x1 matrix
    windVar = (r_bar^2)/3; % 1x1 matrix
    driftVar = 0; % 1x1 matrix
    
    % Initialization
    % x = [px,py,phi,sx,sy,rho,b]
    % v = [vd,vr,vrho,vb]
    % z = [za,zb,zc,zy,zn]
    % w = [wa,wb,wc,wg,wn]
    
    % estimator variance init (initial posterior variance)
    estState.Pm = diag([posVar(1), posVar(2), oriVar, linVelVar(1), linVelVar(2), windVar,  driftVar]);
    % estimator state
    estState.xm = [posEst, oriEst, linVelEst, windEst, driftEst];
    % time of last update
    estState.tm = tm;
end

%% Estimator iteration.

% x = [px,py,phi,sx,sy,rho,b]
% v = [vd,vr,vrho,vb]
% z = [za,zb,zc,zy,zn]
% w = [wa,wb,wc,wg,wn]

% load params
Cdh = estConst.dragCoefficientHydr;
Cda = estConst.dragCoefficientAir;
Cr = estConst.rudderCoefficient;
Cw = estConst.windVel;

Qd = estConst.DragNoise; Qr = estConst.RudderNoise;
Qrho = estConst.WindAngleNoise; Qb = estConst.GyroDriftNoise;
sa = estConst.DistNoiseA; sb = estConst.DistNoiseB;
sc = estConst.DistNoiseC; sg = estConst.GyroNoise;
sn = estConst.CompassNoise;

% get time since last estimator update
dt = tm - estState.tm;
estState.tm = tm; % update measurement update time

% prior update:
% solve y_dot = q(y,0,t) for t = [(k-1)T,kT], y((k-1)T) = xm(k-1)
% solve P_dot = AP + PA' + LQL' for t = [(k-1)T,kT], P((k-1)T) = Pm(k-1)

t_span = [tm-dt tm]; % for t = [(k-1)T,kT]

y0 = estState.xm; % y((k-1)T) = xm(k-1)
P0 = estState.Pm; % P((k-1)T) = Pm(k-1)

% x = [px,py,phi,sx,sy,rho,b]
% v = [vd,vr,vrho,vb]

SQT = @(y) sqrt((y(4) - Cw * cos(y(6)))^2 + (y(5) - Cw * sin(y(6)))^2);

ydot = @(y) [y(4);
             y(5);
             Cr * actuate(2);
             cos(y(3))*(tanh(actuate(1)) - Cdh*( x(4)^2 + x(5)^2)) - Cda*(y(4)-Cw*cos(y(6)))*SQT(y);
             sin(y(3))*(tanh(actuate(1)) - Cdh*( x(4)^2 + x(5)^2)) - Cda*(y(5)-Cw*sin(y(6)))*SQT(y);
             0;
             0];

% A = dq/dx | x=y, v=0
% L = dq/dv | x=y, v=0
A43 = @(y) -sin(y(3))*(tanh(actuate(1)) - Cdh*(y(4)^2+y(5)^2));
A44 = @(y) -cos(y(3))*Cdh*2*y(4) - Cda*SQT(y) - Cda*(y(4)-Cw*cos(y(6)))^2/SQT(y);
A45 = @(y) -cos(y(3))*Cdh*2*y(5) - Cda*(y(4)-Cw*cos(y(6)))*(y(5)-Cw*sin(y(6)))/SQT(y);
A46 = @(y) -Cda*Cw*sin(y(6))*SQT(y) ...
           -Cda*((y(4)-Cw*cos(y(6)))*(Cw*sin(y(6))*(y(4)-Cw*cos(y(6))) - Cw*cos(y(6))*(y(5)-Cw*sin(y(6)))))/SQT(y);
A53 = @(y) cos(y(3))*(tanh(actuate(1)) - Cdh*(y(4)^2 + y(5)^2));
A54 = @(y) -sin(y(3))*Cdh*2*y(4) - Cda*(y(5)-Cw*sin(y(6)))*(y(4)-Cw*cos(y(6)))/SQT(y);
A55 = @(y) -sin(y(3))*Cdh*2*y(5) - Cda*SQT - Cda*(y(5)-Cw*sin(y(6)))^2/SQT(y);
A56 = @(y) Cda*Cw*cos(y(6))*SQT(y) ...
           -Cda*((y(5)-Cw*sin(y(6)))*(Cw*sin(y(6))*(y(4)-Cw*cos(y(6))) - Cw*cos(y(6))*(y(5)-Cw*sin(y(6)))))/SQT(y);
       
A = @(y) [0, 0, 0,      1,      0,      0,      0;
          0, 0, 0,      0,      1,      0,      0;
          0, 0, 0,      0,      0,      0,      0;
          0, 0, A43(y), A44(y), A45(y), A46(y), 0;
          0, 0, A53(y), A54(y), A55(y), A56(y), 0;
          0, 0, 0,      0,      0,      0,      0;
          0, 0, 0,      0,      0,      0,      0;];
      
L = @(y) [0,                              0,             0, 0;
          0,                              0,             0, 0;
          0,                              Cr*actuate(2), 0, 0;
          -cos(y(3))*Cdh*(y(4)^2+y(5)^2), 0,             0, 0;
          -sin(y(3))*Cdh*(y(4)^2+y(5)^2), 0,             0, 0;
          0,                              0,             1, 0;
          0,                              0,             0, 1;];
      
Q = diag([Qd, Qr, Qrho, Qb]);
R = diag([sa, sb, sc, sg, sn]);

Pdot = @(y) A(y)*P(y) + P(y)*A(y)' + L(y)*Q*L(y)';











% [t,y] = ode45(@(t,y) odefcn(t,y,actuate,estConst), t_span, y0);
% xp = y(tm);



% measurement update

% Get resulting estimates and variances
% Output quantities
posEst = estState.xm(1:2);
linVelEst = ...
oriEst = ...
windEst = ...
driftEst = ...

posVar = ...
linVelVar = ...
oriVar = ...
windVar = ...
driftVar = ...

end