function [postParticles] = Estimator(prevPostParticles, sens, act, estConst, km)
% The estimator function. The function will be called in two different
% modes: If km==0, the estimator is initialized. If km > 0, the
% estimator does an iteration for a single sample time interval using the 
% previous posterior particles passed to the estimator in
% prevPostParticles and the sensor measurement and control inputs.
%
% Inputs:
%   prevPostParticles   previous posterior particles at time step k-1
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles) 
%                       corresponding to: 
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]
%                       .kappa: wall offset [m]
%                           
%   sens                Sensor measurement z(k), scalar
%
%   act                 Control inputs u(k-1), [1x2]-vector
%                       act(1): u_f, forward control input
%                       act(2): u_phi, angular control input
%
%   estConst            estimator constants (as in EstimatorConst.m)
%
%   km                  time index k, scalar
%                       corresponds to continous time t = k*Ts
%                       If km==0 initialization, otherwise estimator
%                       iteration step.
%
% Outputs:
%   postParticles       Posterior particles at time step k
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles) 
%                       corresponding to: 
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]
%                       .kappa: wall offset [m]
%
%
% Class:
% Recursive Estimation
% Spring 2021
% Programming Exercise 2
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Matthias Hofer, Carlo Sferrazza
% hofermat@ethz.ch
% csferrazza@ethz.ch

% Set number of particles:
N = 1000; % N_particles, obviously, you will need more particles than 10.

%% Mode 1: Initialization
if (km == 0)
    % Do the initialization of your estimator here!
    % (x,y) uniformly distributed in circle R. 
    % We can prove that (r^2, theta) is also uniformly distributed
    % see blog here: https://www.cnblogs.com/tenosdoit/p/4025221.html
    
    r = estConst.d * sqrt(rand(1,N));
    theta = 2*pi * rand(1,N);
    x = r .* cos(theta);
    y = r .* sin(theta);
    
    xA = estConst.pA(1); yA = estConst.pA(2);
    xB = estConst.pB(1); yB = estConst.pB(2);
    phi0 = estConst.phi_0;
    L = estConst.l;
    
    pos = randi([1,2],1,N) - 1; % 0 for A and 1 for B
    
    postParticles.x_r = x + (1-pos) * xA + pos * xB; % 1xN_particles matrix
    postParticles.y_r = y + (1-pos) * yA + pos * yB; % 1xN_particles matrix
    postParticles.phi = 2*phi0 * rand(1,N) - phi0; % 1xN_particles matrix
    postParticles.kappa = 2*L * rand(1,N) - L; % 1xN_particles matrix
    
    % and leave the function
    return;
end % end init

%% Mode 2: Estimator iteration.
% If km > 0, we perform a regular update of the estimator.
% Implement your estimator here!

% Prior Update:
% x^n_p(k) = q_k?1(x^n_m(k-1),v^n(k-1))
% define v^n(k-1)
v_f   = estConst.sigma_f * (rand(1,N) - 0.5);
v_phi = estConst.sigma_phi * ( rand(1,N) - 0.5);
% define x^n_p(k)
postParticles.x_r = prevPostParticles.x_r + (act(1) + v_f) .* cos(prevPostParticles.phi);
postParticles.y_r = prevPostParticles.y_r + (act(1) + v_f) .* sin(prevPostParticles.phi);
postParticles.phi = prevPostParticles.phi + act(2) + v_phi;
postParticles.kappa = prevPostParticles.kappa;

% Posterior Update:
% scale: --->  beta_n = alpha*p(z|x), sum(beta_n) = 1
% step1: get z by calling get_distance






postParticles.x_r = ...
postParticles.y_r = ...
postParticles.phi = ...
postParticles.kappa = ...

end % end estimator