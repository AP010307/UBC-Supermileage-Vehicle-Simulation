function [energy1,energy2,energy3,time_req,v_end,energy_end] = fcvLapMC(maxV, initialV, distToBrake, initialOutput, numTrials)
% fcvLapMC performs one lap of the FCV monte carlo simulation
%
% Input Parameters:
% maxV: array of normally distributed values for the maximum velocity for the motor to keep running
% initialV: array of values for the initial velocity of the vehicle at the beginning of lap 
% distToBrake: array of normally distributed values for the selected braking distance for the lap
% initialOutput: array of normally distributed values for the intial motor output percent at the beginning
% numTrials: number of trials to run for the monte carlo simulation
%
% Output Parameters
% energy1: total motor input calculated from torque values in kJ
% energy2: total motor input mileage calculated from interpolated motor power in kJ
% energy3: total vehicle mileage calculated from actual motor power values in kJ 
% time_req: total time to complete lap, in seconds
% v_end: the velocity at the end of lap
% energy_end: the vehicle kinetic energy at the end of lap
%
% Read SEMA Track and FCV Motor Data
x = readmatrix('VehicleEfficiencyCalculatorData.xlsx','Sheet','Sonoma_Data','Range','A2:A290');
elev = readmatrix('VehicleEfficiencyCalculatorData.xlsx','Sheet','Sonoma_Data','Range','B2:B290');
torque_data = readmatrix('VehicleEfficiencyCalculatorData.xlsx','Sheet','FCV_Motor','Range','D26:D45'); 
wattsOut_data = readmatrix('VehicleEfficiencyCalculatorData.xlsx','Sheet','FCV_Motor','Range','E26:E45');

% ARRAY SETUP
dist = zeros(1,length(x));
incline_angle = zeros(1,length(x));
incline_angle_sum = zeros(1,length(x));
roll_resis = zeros(1,length(x));
roll_resis_sum = zeros(1,length(x));
pot_energy = zeros(1,length(x));
pot_energy_sum = zeros(1,length(x));
time_req = zeros(1,numTrials);
v_end = zeros(1,numTrials);
energy_end = zeros(1,numTrials);
energy1 = zeros(1,numTrials);
energy2 = zeros(1,numTrials);
energy3 = zeros(1,numTrials);

% Load FCV Vehicle and Motor Parameters
P = load("FCV_Parameters.mat");

% VEHICLE PARAMETERS
m = P.vehicle_weight + P.driver_weight/2.2; % in kg
drag_area = P.frontal_area*P.drag_coeff*P.aero_modif;
rollresis_coeff = P.rollresis_coeff;
wheel_rad = P.wheel_rad; 
density = 1.225;
g = 9.81;

% MOTOR PARAMETERS
supplyVoltage = P.supplyVoltage; 
stallTorque = P.stallTorque; 
noloadSpeed = P.noloadSpeed; 
speedTorqueSlope = P.speedTorqueSlope; 
dampingFactor = P.dampingFactor; 
staticFriction = P.staticFriction; 
Kt = P.Kt; 
windingResis = P.windingResis; 
stallCurrent = supplyVoltage/windingResis; 
gear_ratio = P.gear_ratio;
drivetrain_eff = P.drivetrain_eff;

% Power Output interpolation using polynomial fit
power_polyfit = polyfit(torque_data, wattsOut_data, 2);

for i = 1:length(x)
    dist(i) = sqrt((x(i).^2 + elev(i).^2));
end    

for i = 2:length(x)
    % Calculate incline angle, rolling resistance, potential energy
    incline_angle(i) = rad2deg(atan((elev(i)-elev(i-1))/(dist(i)-dist(i-1))));
    incline_angle_sum(i) = incline_angle_sum(i-1)+incline_angle(i);
    roll_resis(i) = rollresis_coeff*m*g*cosd(incline_angle(i))*(dist(i)-dist(i-1));
    roll_resis_sum(i) = roll_resis_sum(i -1)+roll_resis(i);
    pot_energy(i) = m*g*sind(incline_angle(i))*(dist(i)-dist(i-1));
    pot_energy_sum(i) = pot_energy_sum(i -1)+pot_energy(i);
end

parfor r = 1:numTrials
    velocity = zeros(1,length(x));
    velocity(1) = initialV(r);
    w = zeros(1,length(x)); 
    torque = zeros(1,length(x)); 
    kin_energy = zeros(1,length(x));
    kin_energy(1) = 0.5*m*velocity(1)^2;
    time = zeros(1,length(x));
    aero_drag = zeros(1,length(x));
    aero_drag_sum = zeros(1,length(x));
    motor_output = zeros(1,length(x));
    motorEnergy1 = zeros(1,length(x));
    interp_powerOut = zeros(1,length(x));
    interp_powerIn = zeros(1,length(x));
    motorEnergy2 = zeros(1,length(x));
    motorEnergy2_sum = zeros(1,length(x));
    actual_powerOut = zeros(1,length(x));
    actual_powerIn = zeros(1,length(x));
    motorEnergy3 = zeros(1,length(x));
    motorEnergy3_sum = zeros(1,length(x));
    motorLosses = zeros(1,length(x));
    elecLosses = zeros(1,length(x));
    motor_eff = ones(1,length(x));
    minEnergy = zeros(1,length(x));
    timeLeft = zeros(1,length(x));
    minVel = zeros(1,length(x));
    maxTime = 222.85;
    powerInt = zeros(1,length(x));
    
    % PI Velocity Controls
    u = zeros(1,length(x)+1); % motor output percent
    u(2) = initialOutput(r);
    Kp = maxV(r)/1; % process gain
    Kc = 1/Kp ; % controller gain
    TauI= 3; % integral time constant
    error = zeros(1,length(x));
    intError = zeros(1,length(x));

    for i = 2:length(x)
        % Calculate aerodynamic drag
        aero_drag(i) = 0.5*density*drag_area*velocity(i-1)^2*(dist(i)-dist(i-1));  
        aero_drag_sum(i) = aero_drag_sum(i-1)+aero_drag(i);
        
        if dist(i)>distToBrake(r) && (kin_energy(i-1)>minEnergy(i-1) && velocity(i-1)>minVel(i-1))
            motor_output(i) = 0;
            motorEnergy1(i) = motorEnergy1(i-1);
            motorEnergy2_sum(i) = motorEnergy2_sum(i-1);
            motorEnergy3_sum(i) = motorEnergy3_sum(i-1);
            
            kin_energy(i) = kin_energy(i-1) + motor_output(i) - aero_drag(i)-roll_resis(i)-pot_energy(i); % overall vehicle energy
            minEnergy(i) = sum(roll_resis(i+1:end)) + sum(pot_energy(i+1:end));
            velocity(i) = sqrt(2*kin_energy(i)/m); % vehicle linear velocity
            w(i) = velocity(i)/wheel_rad * gear_ratio; % driveshaft angular velocity
            time(i) = time(i-1) + (dist(i)-dist(i-1))/velocity(i); % time travelled
            timeLeft(i) = maxTime - time(i);
            minVel(i) = 1.5*(x(end)-dist(i)) / timeLeft(i); % safety factor of 1.5
            
        elseif u(i) == 0 % when motor output % is zero
            motorEnergy1(i) = motorEnergy1(i-1);
            motorEnergy2_sum(i) = motorEnergy2_sum(i-1);
            motorEnergy3_sum(i) = motorEnergy3_sum(i-1);
        
            kin_energy(i) = kin_energy(i-1) + motor_output(i) - aero_drag(i)-roll_resis(i)-pot_energy(i); % overall vehicle energy
            minEnergy(i) = sum(roll_resis(i+1:end)) + sum(pot_energy(i+1:end));
            velocity(i) = sqrt(2*kin_energy(i)/m); % vehicle linear velocity
            w(i) = velocity(i)/wheel_rad * gear_ratio; % driveshaft angular velocity
            time(i) = time(i-1) + (dist(i)-dist(i-1))/velocity(i); % time travelled
            timeLeft(i) = maxTime - time(i);
            minVel(i) = 1.5*(x(end)-dist(i)) / timeLeft(i); % safety factor of 1.5
            
            error(i) = maxV(r)-velocity(i); % calculate velocity error value
            intError(i) = intError(i-1) + error(i)*(time(i)-time(i-1)); % calculate integral term
            u(i+1) = Kc*error(i) +Kc/TauI*intError(i); % determine control output (motor output%)
            
            if u(i+1)>1
                u(i+1)=1;
            end
            if u(i+1)<0
                u(i+1)=0;
                intError(i) = 0; % reset integral accumulator
            end

        else 
            torque(i) = (noloadSpeed - w(i-1)) / speedTorqueSlope; 
            if torque(i) > stallTorque
                torque(i) = stallTorque;
            end
            % calculate mechanical and electrical losses
            motorLosses(i) = w(i-1)*(dampingFactor*w(i-1) + staticFriction); 
            elecLosses(i) = (u(i)*torque(i)/Kt)^2 * windingResis; 
            motor_eff(i) = u(i)*torque(i)*w(i-1) / (u(i)*torque(i)*w(i-1)+motorLosses(i)+elecLosses(i));
            motor_eff(2) = 0.5;
            
            motor_output(i) = u(i)*torque(i)*gear_ratio/wheel_rad * (dist(i)-dist(i-1)); % calculate motor energy output using torque values
            motorEnergy1(i) = motorEnergy1(i-1) + (motor_output(i)/(drivetrain_eff*motor_eff(i))); % total motor energy input from torque
            
            kin_energy(i) = kin_energy(i-1) + motor_output(i) - aero_drag(i)-roll_resis(i)-pot_energy(i); % overall vehicle energy
            minEnergy(i) = sum(roll_resis(i+1:end)) + sum(pot_energy(i+1:end));
            
            velocity(i) = sqrt(2*kin_energy(i)/m); % vehicle linear velocity 
            w(i) = velocity(i)/wheel_rad * gear_ratio; % driveshaft angular velocity
            time(i) = time(i-1) + (dist(i)-dist(i-1))/velocity(i); % time travelled
            timeLeft(i) = maxTime - time(i);
            minVel(i) = 1.5*(x(end)-dist(i)) / timeLeft(i); % safety factor of 1.5
            
            % Velocity control
            error(i) = maxV(r)-velocity(i); % calculate velocity error value
            intError(i) = intError(i-1) + error(i)*(time(i)-time(i-1)); % calculate integral/accumulator term
            u(i+1) = Kc*error(i) +Kc/TauI*intError(i) ; % determine velocity control output (motor output%) 
            
            if u(i+1)>=1
                u(i+1)=1;
            end
            if u(i+1)<0
                u(i+1)=0;
               intError(i) = 0; % Reset accumulator 
            end
            
            % Interpolating power out data to find motor energy
            powerInt(i) = power_polyfit(1) *torque(i)^2 + power_polyfit(2)*torque(i) + power_polyfit(3);
            interp_powerOut(i) = u(i)*powerInt(i);
            interp_powerIn(i) = interp_powerOut(i)/(drivetrain_eff*motor_eff(i));
            interp_powerIn(2) = u(2)*12*stallCurrent; % Initial power input during startup
            motorEnergy2(i) = interp_powerIn(i) * (time(i)-time(i-1)); % motor input using interpolated power*time
            motorEnergy2_sum(i) = motorEnergy2_sum(i-1) + motorEnergy2(i);
            
            % Finding power using torque * velocity
            actual_powerOut(i) = u(i)*torque(i) * w(i-1); % calculate power output using torque*angular velocity
            actual_powerIn(i) = actual_powerOut(i) / (drivetrain_eff*motor_eff(i));
            actual_powerIn(2) = u(2)*12*stallCurrent; % Initial power input during startup
            motorEnergy3(i) = actual_powerIn(i) * (time(i)-time(i-1)); % motor input using torque*angular velocity
            motorEnergy3_sum(i) = motorEnergy3_sum(i-1) + motorEnergy3(i);
        end
    end

    time_req(r) = time(end);
    v_end(r) = velocity(end);
    energy_end(r) = kin_energy(end) % Energy wasted at the end
    energy1(r) = motorEnergy1(end)/1000; % Motor Input Energy in kJ (from Force*Dist)
    energy2(r) = sum(motorEnergy2)/1000; % Motor Input Energy in kJ (from intepolated power*time)
    energy3(r) = sum(motorEnergy3)/1000; % Motor Input Energy in kJ (from actual power*time)
end
end