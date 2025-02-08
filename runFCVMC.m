function [T] = runFCVMC(maxV_mean, maxV_sd, brakeDist_mean, brakeDist_sd, initialOutput_mean, initialOutput_sd,numTrials,generatePlots)
% runFCVMC performs the full run of the monte carlo simulation for the Fuel Cell Vehicle

% Input Parameters:
% maxV_mean, maxV_sd: the mean and standard deviation of the maximum velocity for the motor to keep running
% brakeDist_mean, brakeDist_sd: the mean and standard deviation of the selected braking distance to stop vehicle
% initialOutput_mean, initialOutput_sd: the mean and standard deviation of the intial motor output percent to start vehicle
% numTrials: the number of trials for the monte carlo simulation
% generatePlots: select whether to display histogram plots (1 = yes, 0 = no)
%
% Output Parameters
% T: Table of values with monte carlo simulation results (Time, Mileage, Braking Energy)
%

totalLaps = 7;
distToBrake_mean = 1440-brakeDist_mean;

maxV = zeros(1,numTrials);
initialV = zeros(1,numTrials);
distToBrake = zeros(1,numTrials);
initialOutput = zeros(1,numTrials);
brakeDist = zeros(1,numTrials);

totalEnergy1 = zeros(1,numTrials);
totalEnergy2 = zeros(1,numTrials);
totalEnergy3 = zeros(1,numTrials);
totalTime = zeros(1,numTrials);
brakeEnergy = zeros(1,numTrials);
mileage1 = zeros(1,numTrials);
mileage2 = zeros(1,numTrials);
mileage3 = zeros(1,numTrials);

% generate normal distribution of input parameters
for i = 1:length(maxV)
    maxV(i) = normrnd(maxV_mean, maxV_sd);
    distToBrake(i) = normrnd(distToBrake_mean, brakeDist_sd);
    initialOutput(i) = normrnd(initialOutput_mean, initialOutput_sd);
end

tic;
for lapNumber = 1:totalLaps
    [energy1,energy2,energy3,time,v_end,energy_end] = fcvLapMC(maxV, initialV, brakeDist, initialOutput, numTrials);
        
    for j = 1:numTrials
        % sum up energy and time values from previous laps for each trial
        totalEnergy1(j) = totalEnergy1(j) + energy1(j);
        totalEnergy2(j) = totalEnergy2(j) + energy2(j);
        totalEnergy3(j) = totalEnergy3(j) + energy3(j);
        totalTime(j) = totalTime(j) + time(j);

        % set initial velocity for subsequent laps
        initialV(j) = v_end(j);
        % set motor output percent to 0 if there is already kinetic energy
        if initialV(j) > 0
            initialOutput(j) = 0;
        end
        % set the braking distance for the final lap (#7)
        if lapNumber == totalLaps-1
            brakeDist(j) = distToBrake(j);
        end
    end
end
toc;

for j = 1:numTrials
    % calculate mileage values in m/kWh and total time in minutes for ech
    % trial
    mileage1(j) = ((1440*lapNumber)/1609) / (totalEnergy1(j)/3600);
    mileage2(j) = ((1440*lapNumber)/1609) / (totalEnergy2(j)/3600);
    mileage3(j) = ((1440*lapNumber)/1609) / (totalEnergy3(j)/3600);
    totalTime(j) = totalTime(j)/60;
    brakeEnergy(j) = energy_end(j);
        
end
% SIMULATION RESULT
T = table([mean(totalTime);median(totalTime);std(totalTime);min(totalTime);max(totalTime)], ...
    [mean(mileage1); median(mileage1); std(mileage1);min(mileage1);max(mileage1)], ...
    [mean(mileage2); median(mileage2); std(mileage2); min(mileage2); max(mileage2)], ...
    [mean(mileage3); median(mileage3); std(mileage3); min(mileage3); max(mileage3)], ...
    [mean(brakeEnergy); median(brakeEnergy); std(brakeEnergy); min(brakeEnergy); max(brakeEnergy)], ...
    'VariableNames',{'Total Time [min]','Mileage1 [m/kWh]','Mileage2 [m/kWh]','Mileage3 [m/kWh]','Total Braking Energy [J]'},'RowNames',{'Mean','Median','Std Dev','Min','Max'});

if generatePlots == 1
    figure;
    h1 = histogram(totalTime);
    h1.BinEdges = [floor(min(totalTime/10))*10:.5:ceil(max(totalTime/10))*10];
    xlabel('Total Time [min]');
    title('Total Time');
    grid on;
    
    figure;
    h4 = histogram(time);
    h4.BinEdges = [floor(min(time/10))*10:5:ceil(max(time/10))*10];
    xlabel('Lap Time [s]');
    title('Lap Time');
    grid on;
    
    figure;
    h2 = histogram(mileage1);
    h2.BinEdges = [floor(min([mileage1 mileage2]/10))*10:5:ceil(max([mileage1 mileage2]/10))*10];
    hold on;
    h3 = histogram(mileage2);
    h3.BinEdges = [floor(min([mileage1 mileage2]/10))*10:5:ceil(max([mileage1 mileage2]/10))*10];
    legend('Mileage (torque)', 'Mileage (power interpolation)', 'location','northoutside');
    xlabel('Mileage [m/kWh]');
    title('Mileage');
    grid on;
end
end
