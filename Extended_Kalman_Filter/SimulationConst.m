function const = SimulationConst()
% 
% Define the constants used in the simulation.  These constants are not 
% accessible to the estimator.
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

%% Boat dynamics constants
const.dragCoefficientHydr = 0.1; % C_d,h
const.dragCoefficientAir = 0.06;  % C_d,a
const.rudderCoefficient = 2; % C_r
const.windVel = 0.75;  % C_w

% control input constraints
const.MaxThrust = 0.1;
const.MaxRudderCommand = 0.03;

% Minimum time for a segment, in seconds.
% Should be multiple of sampling time.
const.minSegTime = 0.5;

% Maximum time for a segment, in seconds.
% Should be multiple of sampling time.
const.maxSegTime = 2;

%% Radio measurement constants
const.pos_radioA = [-1000 1000]; %[x_a,y_a]
const.pos_radioB = [2000 0]; %[x_b,y_b]
const.pos_radioC = [0 2000]; %[x_c,y_c]

%% Noise properties
% process noise
% defined by its power spectral density Q_
const.DragNoise = 0.1; % Q_d
const.RudderNoise = 0.01; % Q_r
const.WindAngleNoise = 0.01; % Q_rho
const.GyroDriftNoise = 0.01; % Q_b

% measurement noise
% normally distributed with zero mean and variance \sigma_
const.DistNoiseA = 20; % const.DistNoiseA = \sigma_a^2
const.DistNoiseB = 20; % const.DistNoiseB = \sigma_b^2
const.DistNoiseC = 5; % const.DistNoiseC = \sigma_c^2

const.GyroNoise = 0.01; % \sigma_g^2

const.CompassNoise = 0.5; %\sigma_n^2

%% Starting point
% The boat start with equal probility in a cirlce of radius R0 around the 
% origin
const.StartRadiusBound = 10; % R_0

% The initial orientation is uniformly distributed in the range
% [-\bar{\phi}, \bar{\phi}] 
const.RotationStartBound = pi/8; % \bar{\phi}

% The initial wind direction is uniformly distributed in the range
% [-\bar{\rho}, \bar{\rho}] 
const.WindAngleStartBound = pi/16; % \bar{\rho}

% The initial gyro drift is exactly 0
const.GyroDriftStartBound = 0;

%% Times
% Number of samples of the simulation, the total duration of the simulation 
% is then const.N*const.sampleContinuous in seconds.
const.N = 500;

% The sample time for the continuous dynamics, in seconds.
const.sampleContinuous = 0.1;

% The min sample time for distance measurement C, in seconds.
const.sampleDistanceCmin = 0.5;

% The max sample time for distance measurement C, in seconds.
const.sampleDistanceCmax = 2;