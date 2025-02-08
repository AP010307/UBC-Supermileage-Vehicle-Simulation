function [mileage1, mileage2, mileage3, totalTime, brakeEnergy] = runFCV(max_velocity, initial_velocity, brake_distance, initialOutput, plot_fig)
% runFCV performs the full run of the simulation for the Fuel Cell Vehicle
% and returns the mileage, total time, and braking energy
%
% Input Parameters:
% max_velocity: the maximum velocity for the motor to keep running
% inital velocity: the initial velocity at the beginning of lap 1
% brake_distance: the selected braking distance for the last lap
% initialOutput: the intial motor output percent to start vehicle
% plot_fig: select whether to display figures (will only plot figures 
% for first, second, and last lap) (1 = yes, 0 = no)
%
% Output Parameters
% mileage1: total vehicle mileage calculated from torque values in m/kWh
% mileage2: total vehicle mileage calculated from interpolated motor power in m/kWh
% mileage3: total vehicle mileage calculated from actual motor power values in m/kWh 
% totalTime: total time to complete track, in minutes
% brakeEnergy: kinetic energy wasted due to braking in J (final lap)
%

totalLaps = 7;
totalEnergy1 = 0;
totalEnergy2 = 0;
totalEnergy3 = 0;
totalTime = 0;
brakeDist = 0;
plotFigures = plot_fig;
    
for lapNumber = 1:totalLaps
    [energy1,energy2,energy3,time,v_end,energy_end] = fcvLap(max_velocity, initial_velocity,brakeDist,initialOutput,plotFigures);
    % set initial velocity for subsequent laps
    initial_velocity = v_end;
    % set motor output percent to 0 if there is already kinetic energy
    if initial_velocity > 0
        initialOutput = 0;
    end 
    % sum up energy and time values from previous laps
    totalEnergy1 = totalEnergy1 + energy1;
    totalEnergy2 = totalEnergy2 + energy2;
    totalEnergy3 = totalEnergy3 + energy3;
    totalTime = totalTime + time;
    % set the braking distance for the final lap (#7)
    if lapNumber == totalLaps-1
        brakeDist = brake_distance;
    end
    % if user selected to plot figures, will plot figures for the first, second, and last lap
    if plot_fig == 1
        if lapNumber > 1 && lapNumber < 6
            plotFigures = 0;
        else
            plotFigures = 1;
        end
    end
end
% calculate mileage values in m/kWh and total time in minutes
mileage1 = ((1440*lapNumber)/1609) / (totalEnergy1/3600);
mileage2 = ((1440*lapNumber)/1609) / (totalEnergy2/3600);
mileage3 = ((1440*lapNumber)/1609) / (totalEnergy3/3600);
totalTime = totalTime/60;
brakeEnergy = energy_end; 
end
