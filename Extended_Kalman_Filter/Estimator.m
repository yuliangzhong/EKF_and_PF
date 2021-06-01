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
    % x = [px,py,sx,sy,phi,rho,b]
    % v = [vd,vr,vrho,vb]
    % z = [za,zb,zc,zy,zn]
    % w = [wa,wb,wc,wg,wn]
    
    % estimator variance init (initial posterior variance)
    estState.Pm = diag([posVar(1), posVar(2), linVelVar(1), linVelVar(2), oriVar, windVar,  driftVar]);
    % estimator state
    estState.xm = [posEst, linVelEst, oriEst, windEst, driftEst];
    % time of last update
    estState.tm = tm;
    return; % leave function
end

%% Estimator iteration.

% x = [px,py,sx,sy,phi,rho,b]
% v = [vd,vr,vrho,vb]
% z = [za,zb,zc,zy,zn]
% w = [wa,wb,wc,wg,wn]

% get time since last estimator update
dt = tm - estState.tm;
estState.tm = tm; % update measurement update time

% Step1: Prior Update:

% solve x_dot = q(x,0,t) for t = [(k-1)T,kT], x((k-1)T) = xm(k-1)

t_span = linspace(tm-dt,tm,100); % for t = [(k-1)T,kT]
x0 = estState.xm; % x((k-1)T) = xm(k-1) row vector
[tt,xx] = ode45(@(t,x) xODE(t,x,actuate,estConst), t_span, x0');
xp = xx(end,:); % row vector

% solve P_dot = AP + PA' + LQL' for t = [(k-1)T,kT], P((k-1)T) = Pm(k-1)

P0 = estState.Pm; % P((k-1)T) = Pm(k-1)
[TT,PP] = ode45(@(t,P) PODE(t,P,tt,xx,actuate,estConst), t_span, P0(:));
PP_end = PP(end,:);
Pp = reshape(PP_end,[7,7]);

% Step2: Posterior Update / Measurement Update
% given xp,Pp

% load param
xa = estConst.pos_radioA(1); ya = estConst.pos_radioA(2);
xb = estConst.pos_radioB(1); yb = estConst.pos_radioB(2);
xc = estConst.pos_radioC(1); yc = estConst.pos_radioC(2);

H = zeros(5,7);
H(1,1) = (xp(1)-xa)/sqrt((xp(1)-xa)^2 + (xp(2)-ya)^2);
H(1,2) = (xp(2)-ya)/sqrt((xp(1)-xa)^2 + (xp(2)-ya)^2);
H(2,1) = (xp(1)-xb)/sqrt((xp(1)-xb)^2 + (xp(2)-yb)^2);
H(2,2) = (xp(2)-yb)/sqrt((xp(1)-xb)^2 + (xp(2)-yb)^2);
H(3,1) = (xp(1)-xc)/sqrt((xp(1)-xc)^2 + (xp(2)-yc)^2);
H(3,2) = (xp(2)-yc)/sqrt((xp(1)-xc)^2 + (xp(2)-yc)^2);
H(4,5) = 1;
H(4,7) = 1;
H(5,5) = 1;

M = eye(5);

R = diag([estConst.DistNoiseA, estConst.DistNoiseB, estConst.DistNoiseC, estConst.GyroNoise, estConst.CompassNoise]);

hk = [sqrt((xp(1)-xa)^2 + (xp(2)-ya)^2);
      sqrt((xp(1)-xb)^2 + (xp(2)-yb)^2);
      sqrt((xp(1)-xc)^2 + (xp(2)-yc)^2);
      xp(5) + xp(7);
      xp(5);];
      
z = sense'; % column vector

if isinf(z(3))    %clear measurement C
    H(3,:) = [];
    M = eye(4);
    R = diag([estConst.DistNoiseA, estConst.DistNoiseB, estConst.GyroNoise, estConst.CompassNoise]);
    hk(3,:) = [];
    z(3,:) = [];
end

% Update
% p_size = size(Pp);
% h_size = size(H);
% m_size = size(M);
% r_size = size(R);
% disp('-----------')
% disp(tm)
% disp(p_size)
% disp(h_size)
% disp(m_size)
% disp(r_size)
K = Pp*H'/(H*Pp*H' + M*R*M');
xm = (xp' + K*(z-hk))'; % row vector
Pm = (eye(7) - K*H)*Pp;

% Step3, almost done!
% Get resulting estimates and variances
% Output quantities
estState.xm = xm;
estState.Pm = Pm;
posEst = xm(1:2);
linVelEst = xm(3:4);
oriEst = xm(5);
windEst = xm(6);
driftEst = xm(7);

posVar = [Pm(1,1), Pm(2,2)];
linVelVar = [Pm(3,3), Pm(4,4)];
oriVar = Pm(5,5);
windVar = Pm(6,6);
driftVar = Pm(7,7);
end

%% ODE functions

function dxdt = xODE(t,x,u,estConst)

    % load params
    Cdh = estConst.dragCoefficientHydr;
    Cda = estConst.dragCoefficientAir;
    Cr = estConst.rudderCoefficient;
    Cw = estConst.windVel;

    SQT = sqrt((x(3) - Cw * cos(x(6)))^2 + (x(4) - Cw * sin(x(6)))^2);
    dxdt = zeros(7,1);
    dxdt(1) = x(3);
    dxdt(2) = x(4);
    dxdt(3) = cos(x(5))*(tanh(u(1)) - Cdh*(x(3)^2 + x(4)^2)) - Cda*(x(3)-Cw*cos(x(6)))*SQT;
    dxdt(4) = sin(x(5))*(tanh(u(1)) - Cdh*(x(3)^2 + x(4)^2)) - Cda*(x(4)-Cw*sin(x(6)))*SQT;
    dxdt(5) = Cr * u(2);
end

function dPdt = PODE(t,P,tt,xx,u,estConst)
    % load params
    Cdh = estConst.dragCoefficientHydr;
    Cda = estConst.dragCoefficientAir;
    Cr = estConst.rudderCoefficient;
    Cw = estConst.windVel;
    
    % get x(t)
    x = interp1(tt,xx,t); % interp_one
    
    % get sqrt
    SQT = sqrt((x(3) - Cw * cos(x(6)))^2 + (x(4) - Cw * sin(x(6)))^2);
    
    % get A(t), A = dq/dx | x=x(t), v=0
    A = zeros(7,7);
    A(1,3) = 1;
    A(2,4) = 1;
    A(3,3) = -cos(x(5))*Cdh*2*x(3) - Cda*SQT - Cda*(x(3)-Cw*cos(x(5)))^2/SQT;
    A(3,4) = -cos(x(5))*Cdh*2*x(4) - Cda*(x(3)-Cw*cos(x(6)))*(x(4)-Cw*sin(x(6)))/SQT;
    A(3,5) = -sin(x(5))*(tanh(u(1)) - Cdh*(x(3)^2+x(4)^2));
    A(3,6) = -Cda*Cw*sin(x(6))*SQT ...
             -Cda*((x(3)-Cw*cos(x(6)))*(Cw*sin(x(6))*(x(3)-Cw*cos(x(6))) - Cw*cos(x(6))*(x(4)-Cw*sin(x(6)))))/SQT;
    A(4,3) = -sin(x(5))*Cdh*2*x(3) - Cda*(x(4)-Cw*sin(x(6)))*(x(3)-Cw*cos(x(6)))/SQT;
    A(4,4) = -sin(x(5))*Cdh*2*x(4) - Cda*SQT - Cda*(x(4)-Cw*sin(x(6)))^2/SQT;
    A(4,5) = cos(x(5))*(tanh(u(1)) - Cdh*(x(3)^2 + x(4)^2));
    A(4,6) = Cda*Cw*cos(x(6))*SQT ...
            -Cda*((x(4)-Cw*sin(x(6)))*(Cw*sin(x(6))*(x(3)-Cw*cos(x(6))) - Cw*cos(x(6))*(x(4)-Cw*sin(x(6)))))/SQT;
    
    % get L(t), L = dq/dv | x=x(t), v=0
    L = zeros(7,4);
    L(5,2) = Cr*u(2);
    L(3,1) = -cos(x(5))*Cdh*(x(3)^2+x(4)^2);
    L(4,1) = -sin(x(5))*Cdh*(x(3)^2+x(4)^2);
    L(5,2) = Cr*u(2);
    L(6,3) = 1;
    L(7,4) = 1;
    
    % get Q
    Q = diag([estConst.DragNoise, estConst.RudderNoise, estConst.WindAngleNoise, estConst.GyroDriftNoise]);
    
    % P_dot
    PP = reshape(P,[7,7]);
    P_dot = A*PP + PP*A'+L*Q*L';
    dPdt = P_dot(:); % reshape
end
