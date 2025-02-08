%% Supermileage Vehicle Simulation
% Documentation: https://docs.google.com/document/d/1g9IYxAuLaUAbttDUqafE3u3H9sW9519EXk1D3lMdIaE/edit?usp=sharing
% Author: Ambrose Lee
% Last Updated: Aug 22, 2020

% Use this file to set vehicle/motor parameters and input parameters
%#ok<*NASGU>
%% Set FCV Vehicle/Motor Parameters
clear;
save FCV_Parameters;

% VEHICLE PARAMETERS
vehicle_weight = 40; % in kg 
driver_weight = 130; % in lbs
frontal_area = 0.3 ; % in m^2 
drag_coeff = 0.13;
aero_modif = 1.5; 
rollresis_coeff = 0.0023; 
wheel_rad = 0.500/2; % in m

% MOTOR PARAMETERS
% values in () are the conversion factors to SI units from the original data
supplyVoltage = 12; % in V
stallTorque = 840 / (141.61); % in oz-in, converted to N-m
noloadSpeed = 4200 * (2*pi/60); % in rpm, converted to in rad/s
speedTorqueSlope = 5 * (2*pi/60*141.61); % in rpm/oz-in, converted to in rad/s / N-m
dampingFactor = 0.12 / 141.61 * 60 / (2*pi*1000); % in oz-in/krpm, converted to in N-m / rad/s
staticFriction = 1.0 / (141.61); % in oz-in, converted to N-m
windingResis = 0.049; % in ohms
Kt = 3.86 / (141.61); % in oz-in/amp, converted to N-m/amp

gear_ratio = 12.25; 
drivetrain_eff = 0.95;

save FCV_Parameters;
%% Set FCV Simulation Input Parameters
clear;
save FCV_Input;

% REGULAR SIMULATION
% Maximum velocity for the motor to keep running (in m/s)
max_velocity = 4.5;
% Braking distance (final lap) (in m)
brake_distance = 1440;
% Initial motor output percent (first lap)
initialOutput = 0.5;
% Choose to display figures (will only display first, second, last lap) (1 = yes, 0 = no)
plot_fig = 1;

% MONTE CARLO SIMULATION
% Set max velocity mean and standard deviation (in m/s)
maxV_mean = 5.5;
maxV_sd = .5;
% Set braking distance mean and standard deviation (in meters)
brakeDist_mean = 1000;
brakeDist_sd = 50;
% Set initial motor output percent mean and standard deviation
initialOutput_mean = 0.5;
initialOutput_sd = 0.00;
% Set number of runs
numTrials = 100;
% Choose to display histograms (1 = yes, 0 = no)
generateHistogram = 1;

save FCV_Input;

%% Set Urban Vehicle/Motor Parameters

clear;
save Urban_Parameters;

% VEHICLE PARAMETERS
vehicle_weight = 82; % in kg
driver_weight = 130; % in lbs
frontal_area = 1.7; % in m^2
drag_coeff = 0.11;
aero_modif = 2;
rollresis_coeff = 0.0013;
wheel_rad = 0.3305; % in m

% MOTOR PARAMETERS
% values in () are the conversion factors to SI units from the original data
supplyVoltage = 12; % in V
stallTorque = 413 / (141.611); % in oz-in, converted to N-m
noloadSpeed = 8128 *(2*pi/60); % in rpm, converted to in rad/s
dampingFactor = 0.95 / 141.61 * 60 / (2*pi*1000); % in oz-in/krpm, converted to in N-m / rad/s
staticFriction = 2.7 / (141.61); % in oz-in, converted to N-m
Kt = 7.99/ (141.61); % in oz-in/amp, converted to N-m/amp
windingResis = 0.053; % in ohms

gear_ratio = 25;
gearbox_eff = 0.98;

save Urban_Parameters;

%% Set Urban Simulation Input Parameters
clear;
save Urban_Input;

% REGULAR SIMULATION
% Maximum velocity for the motor to keep running (in m/s)
max_velocity = 7;
% Braking distance (final lap) (in m)
brake_distance = 235;
% Initial motor output percent (first lap)
initialOutput = 1;
% Choose to display figures (will only display first lap) (1 = yes, 0 = no)
plot_fig = 0;

% MONTE CARLO SIMULATION
% Set max velocity mean and standard deviation (in m/s)
maxV_mean = 7.5;
maxV_sd = 0.5;
% Set braking distance mean and standard deviation (in meters)
brakeDist_mean = 250;
brakeDist_sd = 15;
% Set initial motor output percent mean and standard deviation
initialOutput_mean = 1;
initialOutput_sd = 0.00;
% Set number of runs
numTrials = 100;
% Choose to display histograms (1 = yes, 0 = no)
generateHistogram = 1;

save Urban_Input;
