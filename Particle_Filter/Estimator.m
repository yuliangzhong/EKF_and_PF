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

%% Configuration here

Sampling_globally = 1;
Sampling_locally = 0;
right_equal = 2;
% Set number of particles:
N = 1000; % N_particles 10k->16min;
remedy_method = 5; % Sampling_locally; % Sampling_globally / Sampling_locally / right_equal
global_ratio = 10; % massive sampling ratio if all of particles fail, if sampling globally
local_variance = 0.2/3; % sample locally -> N(0,localvariance^2), s.t. 3sigma<0.2
K = 0.01; % roughening param


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
% x^n_p(k) = q_k-1(x^n_m(k-1),v^n(k-1))
% define v^n(k-1)
v_f   = estConst.sigma_f * (rand(1,N) - 0.5);
v_phi = estConst.sigma_phi * (rand(1,N) - 0.5);
% define x^n_p(k)
postParticles.x_r = prevPostParticles.x_r + (act(1) + v_f) .* cos(prevPostParticles.phi);
postParticles.y_r = prevPostParticles.y_r + (act(1) + v_f) .* sin(prevPostParticles.phi);
postParticles.phi = prevPostParticles.phi + act(2) + v_phi;
postParticles.kappa = prevPostParticles.kappa;

% Posterior Update:
% scale: --->  beta_n = alpha*p(sens|xp), sum(beta_n) = 1
z = zeros(1,N);
p_sens = zeros(1,N);
for k = 1:N
    z(k) = get_distance(postParticles.x_r(k),postParticles.y_r(k),postParticles.phi(k),postParticles.kappa(k),estConst.contour);
    p_sens(k) = pdf(sens - z(k), estConst.epsilon);
end
SUM = sum(p_sens);
% what if sum(p_sens)=0 ?
% that means num of particles is not big enough
% Remedy with more intensive sampling!
% call intensive_sampling()
while SUM==0
    disp("Warning, sum(p_sens) = 0, intensive resampling! At time = "+num2str(km))
    
    if remedy_method == Sampling_globally
        M = global_ratio*N; % intensive sampling number
        % sample kappa: 2 methods
        kappa = 0.2*(2 * estConst.l * rand(1,M) - estConst.l); % 0.2*U[-L,L]
        % uniformly sample x and y in rectangle
        posx = (3-kappa) .* rand(1,M) + kappa; 
        posy = 3 * rand(1,M);                
        % discard samples outside map
        for i=1:M
            if ~Is_inmap(posx(i),posy(i))
                posx(i) = nan;
                posy(i) = nan;
            end
        end
        % resampling result
        index = ~isnan(posx);
        posx = posx(index);
        posy = posy(index);
        kappa = kappa(index);
        len = length(posx);
        phi = 2*pi* rand(1,len)-pi;
        % calculate p(z|x)
        z_tmp = zeros(1,len);
        p_tmp = zeros(1,len);
        for k = 1:len
            z_tmp(k) = get_distance(posx(k),posy(k),phi(k),kappa(k),estConst.contour);
            p_tmp(k) = pdf(sens - z_tmp(k), estConst.epsilon);
        end
        [p_sens,Ind] = maxk(p_tmp,N);
        postParticles.x_r = posx(Ind);
        postParticles.y_r = posy(Ind);
        postParticles.phi = phi(Ind);
        postParticles.kappa = kappa(Ind);
        
        
    elseif remedy_method == Sampling_locally
        postParticles.x_r = postParticles.x_r + local_variance*randn(1,N);
        postParticles.y_r = postParticles.y_r + local_variance*randn(1,N);
        for k = 1:N
            z(k) = get_distance(postParticles.x_r(k),postParticles.y_r(k),postParticles.phi(k),postParticles.kappa(k),estConst.contour);
            p_sens(k) = pdf(sens - z(k), estConst.epsilon);
        end
        
        
    elseif remedy_method == right_equal % not recommended
        p_sens = 1/N * ones(1,N);
    else
        disp("no remedy method, failure!")
        break
    end
    SUM = sum(p_sens);
end
beta_n = p_sens / SUM; % i.e. --> beta_n

% resampling from old particles
% xm = resampling(beta_n);
cumSUM = cumsum(beta_n);
ind = zeros(1,N);
for i = 1:N
    r = rand();
    j = N;
    while r<= cumSUM(j-1)
        j = j-1;
        if j==1
            break;
        end
    end
    ind(i) = j;
end
for k = 1:N
    postParticles.x_r(k) = postParticles.x_r(ind(k));
    postParticles.y_r(k) = postParticles.y_r(ind(k));
    postParticles.phi(k) = postParticles.phi(ind(k));
    postParticles.kappa(k) = postParticles.kappa(ind(k));
end

% roughening: Lecture11 Page23
Ex = max(postParticles.x_r) - min(postParticles.x_r);
Ey = max(postParticles.y_r) - min(postParticles.y_r);
Ep = max(postParticles.phi) - min(postParticles.phi);
Ek = max(postParticles.kappa) - min(postParticles.kappa);

sigma_x = K * Ex * N^(-1/4);
sigma_y = K * Ey * N^(-1/4);
sigma_p = K * Ep * N^(-1/4);
sigma_k = K * Ek * N^(-1/4);
postParticles.x_r = postParticles.x_r + sigma_x*randn(1,N);
postParticles.y_r = postParticles.y_r + sigma_y*randn(1,N);
postParticles.phi = postParticles.phi + sigma_p*randn(1,N);
postParticles.kappa = postParticles.kappa + sigma_k*randn(1,N);

end % end estimator


% how to know whether a ray and a line segment intersect?
% the algorithm is adapted from the website below
% https://stackoverflow.com/questions/14307158/
function distance = get_distance(xr,yr,phi,kappa,contour)
    contour(8,1) = kappa;
    contour(9,1) = kappa;
    contour = [contour; contour(1,:)];
    distance = 99999;
    for seg = 1:10
        % (x1,y1) --(x,y)-- (x2,y2)
        x1 = contour(seg,1); y1 = contour(seg,2);
        x2 = contour(seg+1,1); y2 = contour(seg+1,2);
        % (x,y) = (x1,y1) + u(x2-x1,y2-y1)= (xr,yr) + t(cos(phi),sin(phi))
        A = [x2 - x1, -cos(phi);
             y2 - y1, -sin(phi)];
        b = [xr - x1;
             yr - y1];
        if rank(A)~=rank([A,b]) % no solution
        elseif rank(A) < 2 && dot([x1-xr; y1-yr],[cos(phi);sin(phi)])>0 % inf solutions => collinear
            distance = min([distance,norm([x1-xr; y1-yr]),norm([x2-xr; y2-yr])]);
        elseif rank(A) == 2 % unique solution
            C = A\b;
            if C(1)>=0 && C(1)<=1 && C(2)>=0
                xc = x1 + C(1)*(x2 - x1);
                yc = y1 + C(1)*(y2 - y1);
                distance = min(distance,norm([xc-xr; yc-yr]));
            end
        end
    end
end

function pw = pdf(w,eps)
    w = abs(w);
    if w < 2*eps
        pw = -w/(5*eps^2) + 2/5/eps;
    elseif w<2.5*eps
        pw = 2/(5*eps^2)*(w - 2*eps);
    elseif w<3*eps
        pw = -2/(5*eps^2)*(w - 3*eps);
    else
        pw = 0;
    end
end

function result = Is_inmap(x,y)
    cond1 = x<=0.5 && y<=0.5;
    cond2 = x<=1 && y>=2.5;
    cond3 = x>=1 && x<=1.25 && y>=-x+3.5;
    cond4 = x>=1.25 && x<=2 && y>=x+1;
    cond5 = x>=2 && y>= -x+5;
    cond6 = x>=2.5 && y<= x-1;
    if cond1||cond2||cond3||cond4||cond5||cond6
        result = false;
    else
        result = true;
    end
end
