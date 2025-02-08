function [mileage1, mileage2, mileage3, totalTime, brakeEnergy] = runUrban(max_velocity, brake_distance, initialOutput, plot_fig)
% runUrban performs the full run of the simulation for Urban Concept
% vehicle
%
% Input Parameters:
% max_velocity: the maximum velocity for the motor to keep running
% brake_distance: the selected braking distance to stop vehicle
% initialOutput: the intial motor output percent to start vehicle
% plot_fig: select whether to display figures (1 = yes, 0 = no)
%
% Output Parameters
% mileage1: total vehicle mileage calculated from torque values in m/kWh
% mileage2: total vehicle mileage calculated from interpolated motor power in m/kWh
% mileage3: total vehicle mileage calculated from actual motor power values in m/kWh 
% totalTime: total time to complete track, in minutes
% brakeEnergy: total kinetic energy wasted due to braking
%

totalLaps = 7;
totalEnergy1 = 0;
totalEnergy2 = 0;
totalEnergy3 = 0;
totalTime = 0;
brakeEnergy = 0;
plotFigures = plot_fig;
    
for lapNumber = 1:totalLaps
    [energy1,energy2,energy3,time,energy_end] = urbanLap(max_velocity,brake_distance,initialOutput,plotFigures);
    % sum up energy and time values from previous laps
    totalEnergy1 = totalEnergy1 + energy1;
    totalEnergy2 = totalEnergy2 + energy2;
    totalEnergy3 = totalEnergy3 + energy3;
    totalTime = totalTime + time;
    brakeEnergy = brakeEnergy + energy_end;
    % if user selected to plot figures, will plot figures for the first lap only
    plotFigures = 0;
end
% calculate mileage values in m/kWh and total time in minutes
mileage1 = ((1440*lapNumber)/1609) / (totalEnergy1/3600);
mileage2 = ((1440*lapNumber)/1609) / (totalEnergy2/3600);
mileage3 = ((1440*lapNumber)/1609) / (totalEnergy3/3600);
totalTime = totalTime/60;
end