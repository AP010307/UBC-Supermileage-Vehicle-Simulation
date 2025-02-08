function [energy1,energy2,energy3,time_req,energy_end] = urbanLap(max_velocity, brake_distance, initialOutput, plotFigures)
% urbanLap performs one lap of the Urban Concept vehicle simulation
%
% Input Parameters:
% max_velocity: the maximum velocity for the motor to keep running
% brake_distance: the selected braking distance for the lap
% initialOutput: the intial motor output percent at the beginning
% plotFigures: select whether to display figures (1 = yes, 0 = no)
%
% Output Parameters
% energy1: total motor input calculated from torque values in kJ
% energy2: total motor input mileage calculated from interpolated motor power in kJ
% energy3: total vehicle mileage calculated from actual motor power values in kJ 
% time_req: total time to complete lap, in seconds
% energy_end: the vehicle kinetic energy at the end of lap
%

% Read SEMA Track and Urban Motor Data
x = readmatrix('VehicleEfficiencyCalculatorData.xlsx','Sheet','Sonoma_Data','Range','A2:A290');
elev = readmatrix('VehicleEfficiencyCalculatorData.xlsx','Sheet','Sonoma_Data','Range','B2:B290');
torque_data = readmatrix('VehicleEfficiencyCalculatorData.xlsx','Sheet','Urban_Motor','Range','D23:D34');
rpm_data = readmatrix('VehicleEfficiencyCalculatorData.xlsx','Sheet','Urban_Motor','Range','B23:B34');
wattsOut_data = readmatrix('VehicleEfficiencyCalculatorData.xlsx','Sheet','Urban_Motor','Range','E23:E34');

% Load Urban Vehicle and Motor Parameters
P = load("Urban_Parameters.mat");

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
dampingFactor = P.dampingFactor; 
staticFriction = P.staticFriction; 
Kt = P.Kt; 
windingResis = P.windingResis; 
stallCurrent = supplyVoltage/windingResis; 

gear_ratio = P.gear_ratio;
gearbox_eff = P.gearbox_eff;

% PI Velocity Controller Parameters
Kp = max_velocity/1; % process gain
Kc = 1/Kp ; % controller gain
TauI= 3; % integral time constant
    
% ARRAY SETUP
dist = zeros(1,length(x));
incline_angle = zeros(1,length(x));
incline_angle_sum = zeros(1,length(x));
velocity = zeros(1,length(x));
velocity(1) = 0;
w = zeros(1,length(x));
w(1) = velocity(1)/wheel_rad * gear_ratio;
torque = zeros(1,length(x));
kin_energy = zeros(1,length(x));
kin_energy(1) = 0.5*m*velocity(1)^2;
time = zeros(1,length(x));
aero_drag = zeros(1,length(x));
aero_drag_sum = zeros(1,length(x));
roll_resis = zeros(1,length(x));
roll_resis_sum = zeros(1,length(x));
pot_energy = zeros(1,length(x));
pot_energy_sum = zeros(1,length(x));
motor_output = zeros(1,length(x));
motorEnergy_1 = zeros(1,length(x));
interp_powerOut = zeros(1,length(x));
interp_powerIn = zeros(1,length(x));
motorEnergy2 = zeros(1,length(x));
motorEnergy2_sum = zeros(1,length(x));
actual_powerOut = zeros(1,length(x));
actual_powerIn = zeros(1,length(x));
motorEnergy3 = zeros(1,length(x));
motorEnergy3_sum = zeros(1,length(x));
motor_eff = zeros(1,length(x));
motorLosses = zeros(1,length(x));
elecLosses = zeros(1,length(x));
power_interp = zeros(1,length(x));
distToBrake = x(end) - brake_distance;
u = zeros(1,length(x)+1); % motor output percent
u(2) = initialOutput;
error = zeros(1,length(x));
intError = zeros(1,length(x));
minEnergy = zeros(1,length(x));
timeLeft = zeros(1,length(x));
minVel = zeros(1,length(x));
maxTime = 222.85;

% Motor Torque + Angular Velocity interpolation using polynomial fit
p = polyfit(torque_data, rpm_data, 1); % assume linear relation
% Motor Power Out interpolation
power_polyfit = polyfit(torque_data, wattsOut_data, 2); % assume quadratic relation

% Simulation Start
for i = 1:length(x)
    dist(i) = sqrt((x(i).^2 + elev(i).^2)); % determine actual distance travelled including elevation changes
end
for i = 2:length(x)
    % Calculate incline angle, rolling resistance and potential energy
    incline_angle(i) = rad2deg(atan((elev(i)-elev(i-1))/(dist(i)-dist(i-1))));
    incline_angle_sum(i) = incline_angle_sum(i-1)+incline_angle(i);
    
    roll_resis(i) = rollresis_coeff*m*g*cosd(incline_angle(i))*(dist(i)-dist(i-1));
    roll_resis_sum(i) = roll_resis_sum(i -1)+roll_resis(i);
    
    pot_energy(i) = m*g*sind(incline_angle(i))*(dist(i)-dist(i-1));
    pot_energy_sum(i) = pot_energy_sum(i -1)+pot_energy(i);
end

for i = 2:length(x)
    % Calculate aerodynamic drag
    aero_drag(i) = 0.5*density*drag_area*velocity(i-1)^2*(dist(i)-dist(i-1));
    aero_drag_sum(i) = aero_drag_sum(i-1)+aero_drag(i);

    % check if vehicle is past the braking distance
    if dist(i)>distToBrake && (kin_energy(i-1)>minEnergy(i-1) && velocity(i-1)>minVel(i-1))
        motor_output(i) = 0; % motor off
        motorEnergy_1(i) = motorEnergy_1(i-1);
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
        motorEnergy_1(i) = motorEnergy_1(i-1);
        motorEnergy2_sum(i) = motorEnergy2_sum(i-1);
        motorEnergy3_sum(i) = motorEnergy3_sum(i-1);
        
        kin_energy(i) = kin_energy(i-1) + motor_output(i) - aero_drag(i)-roll_resis(i)-pot_energy(i); % overall vehicle energy
        minEnergy(i) = sum(roll_resis(i+1:end)) + sum(pot_energy(i+1:end));
        velocity(i) = sqrt(2*kin_energy(i)/m); % vehicle linear velocity
        w(i) = velocity(i)/wheel_rad * gear_ratio; % driveshaft angular velocity
        time(i) = time(i-1) + (dist(i)-dist(i-1))/velocity(i); % time travelled
        timeLeft(i) = maxTime - time(i);
        minVel(i) = 1.5*(x(end)-dist(i)) / timeLeft(i); % safety factor of 1.5
        
        error(i) = max_velocity-velocity(i); % calculate velocity error value
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
        torque(i) = fzero(@(x)  p(1) .*x + p(2) - w(i-1), 1); % find torque based on the current motor angular velocity
        if torque(i) > stallTorque
            torque(i) = stallTorque;
        end
        
        % calculate mechanical and electrical losses and motor efficiency
        motorLosses(i) = w(i-1)*(dampingFactor*w(i-1) + staticFriction); 
        elecLosses(i) = (u(i)*torque(i)/Kt)^2 * windingResis; 
        motor_eff(i) = u(i)*torque(i)*w(i-1) / (u(i)*torque(i)*w(i-1)+motorLosses(i)+elecLosses(i));
        motor_eff(2) = 0.5; % assumed motor eff (starting value)
        
        motor_output(i) = u(i)*torque(i)*gear_ratio/wheel_rad * (dist(i)-dist(i-1)); % calculate motor energy output using torque values
        motorEnergy_1(i) = motorEnergy_1(i-1) + (motor_output(i)/(gearbox_eff*motor_eff(i))); % total motor energy input from torque
        
        kin_energy(i) = kin_energy(i-1) + motor_output(i) - aero_drag(i)-roll_resis(i)-pot_energy(i); % overall vehicle energy
        minEnergy(i) = sum(roll_resis(i+1:end)) + sum(pot_energy(i+1:end));
        velocity(i) = sqrt(2*kin_energy(i)/m); % vehicle linear velocity 
        w(i) = velocity(i)/wheel_rad * gear_ratio; % driveshaft angular velocity
        time(i) = time(i-1) + (dist(i)-dist(i-1))/velocity(i); % time travelled
        timeLeft(i) = maxTime - time(i);
        minVel(i) = 1.5*(x(end)-dist(i)) / timeLeft(i); % safety factor of 1.5
        
        % Velocity control
        error(i) = max_velocity-velocity(i); % calculate velocity error value
        intError(i) = intError(i-1) + error(i)*(time(i)-time(i-1)); % calculate integral/accumulator term
        u(i+1) = Kc*error(i) + Kc/TauI*intError(i) ; % determine velocity control output (motor output%) 
        
        if u(i+1)>=1
            u(i+1)=1;
        end
        if u(i+1)<0
            u(i+1)=0;
           intError(i) = 0; % Reset accumulator 
        end
        
        % Interpolating power out data to find motor energy
        power_interp(i) = power_polyfit(1) *torque(i)^2 + power_polyfit(2)*torque(i) + power_polyfit(3);
        interp_powerOut(i) = u(i)*power_interp(i);
        interp_powerIn(i) = interp_powerOut(i)/(gearbox_eff*motor_eff(i));
        %interp_powerIn(2) = u(2)*supplyVoltage*stall_current; % Initial power input during start
        motorEnergy2(i) = interp_powerIn(i) * (time(i)-time(i-1)); % motor input using interpolated power*time
        motorEnergy2_sum(i) = motorEnergy2_sum(i-1) + motorEnergy2(i);
        
        % Finding power using torque * velocity
        actual_powerOut(i) = u(i)*torque(i) * w(i); % calculate power output using torque*angular velocity
        actual_powerIn(i) = actual_powerOut(i) / (gearbox_eff*motor_eff(i));
        %actual_powerIn(2) = u(2)*supplyVoltage*stall_current; % Initial power input during start
        motorEnergy3(i) = actual_powerIn(i) * (time(i)-time(i-1)); % motor input using torque*angular velocity
        motorEnergy3_sum(i) = motorEnergy3_sum(i-1) + motorEnergy3(i);
    end
end

% adjusting PE values to make the minimum PE = zero
minPE = abs(min(pot_energy_sum));
for m = 1:length(x)
    pot_energy_sum(m) = pot_energy_sum(m) + minPE;
end

% SIMULATION RESULTS %
time_req = time(end); % Time required to complete lap
v_end = velocity(end); % Remaining vehicle velocity
energy_end = kin_energy(end); % Excess energy wasted
v_max = max(velocity); % Max vehicle velocity throughout lap
energy1 = motorEnergy_1(end)/1000; % Motor Input Energy in kJ (from Force*Dist)
energy2 = sum(motorEnergy2)/1000; % Motor Input Energy in kJ (from intepolated power*time)
energy3 = sum(motorEnergy3)/1000; % Motor Input Energy in kJ (from actual power*time)

lapMileage1 = (x(end)/1609) / (motorEnergy_1(end)/1000/3600); % Mileage in m/kWh (from Force*Dist)
lapMileage2 = (x(end)/1609) / (sum(motorEnergy2)/1000/3600); % Mileage in m/kWh (from intepolated power*time)
lapMileage3 = (x(end)/1609) / (sum(motorEnergy3)/1000/3600); % Mileage in m/kWh (torque*velocity*time)

if plotFigures == 1
    % Vehicle Energy Profile Figure
    figure;
    Y = [transpose(kin_energy), transpose(aero_drag_sum), transpose(roll_resis_sum), transpose(pot_energy_sum) ];
    area(transpose(dist),Y)
    xlim([0 1500]);  xlabel("Distance from Start [m]"); ylabel("Energy [J]")
    legend('Kinetic Energy','Aerodynamic Drag','Rolling Resistance','Potential Energy','Location', 'southoutside');
    title ('Vehicle Energy Profile')
    
    % Velocity/Elevation/Motor Input Profile
    figure;
    yyaxis left; 
    plot(dist,elev,'b-', dist, velocity, 'r-')
    yline(max_velocity, 'k--');
    xlabel ("Distance from Start [m]"); ylabel ("Elevation [m], Velocity [m/s]"); 
    yyaxis right; 
    plot(dist, motorEnergy_1, 'g-', dist, motorEnergy2_sum, 'k-')
    ylabel("Energy [J]");
    grid on; xlim([0 1500]); legend('Elevation Change from Start','Velocity','Max Velocity Set','Motor Energy (from Torque)', 'Motor Energy (Interpolated Power values)','Location', 'southoutside'); legend('boxoff');
    title ('Velocity/Elevation/Motor Energy Profile');
    
    % Power Consumption and Motor Efficiency
    figure;
    plot(dist, actual_powerOut, dist, actual_powerIn);
    xlabel("Distance [m]")
    ylabel("Power [W]")
    yyaxis right;
    plot(dist, motor_eff, 'k.');
    ylabel("Motor Efficiency")
    legend('Power Out', 'Power In', 'Motor Efficiency', 'Location', 'southoutside');
    title ('Power Consumption and Motor Efficiency');
    
    % Velocity vs Time
    figure;
    plot(time, velocity) 
    xlabel ("Time [s]"); ylabel ("Velocity [m/s]"); 
    title("Velocity vs. Time")
    grid on;
    
    % PI Controller Output
    figure;
    plot(dist,elev,'b-', dist, velocity, 'r-')
    yyaxis right;
    plot(dist, u(2:end), 'k');
    ylabel("Motor Output %");
    legend('Elevation','Velocity','Motor Output %','Location', 'southoutside');
    title ('PI Controller Output (Motor output percent)')
    
    % PI Controller Error     
%     figure;
%     plot(dist, velocity, 'b',dist, error, 'r', dist, intError, 'g')
%     hold on;
%     yline(max_velocity);
%     yline(0);
%     legend('Velocity', 'Velocity Error', 'Integral of Error','Location', 'southoutside')
%     grid on;
%     title ('PI Controller Output (Velocity Error)')
end
end
