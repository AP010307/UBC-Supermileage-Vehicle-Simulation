%% Load Data
track_data = readmatrix('sem_2023_us.csv');
x = track_data(:, 1);
y = track_data(:, 2);
z = track_data(:, 3);

% Check Data
disp('X, Y, Z Data:');
disp([x, y, z]);

%% 2D Track Visualization
figure;
plot(x, y, 'k-', 'LineWidth', 2); % 2D track plot
grid on;
title('2D Track Visualization');
xlabel('X (meters)');
ylabel('Y (meters)');
axis equal; % Maintain aspect ratio

set(gca, 'XDir', 'reverse'); % Reverse the direction of the Y-axis

disp(['X range: ', num2str(min(x)), ' to ', num2str(max(x))]);
disp(['Y range: ', num2str(min(y)), ' to ', num2str(max(y))]);
disp(['Z range: ', num2str(min(z)), ' to ', num2str(max(z))]);


%% Scale Factors
scale_factor = 1000; % Adjust based on your needs
x_scaled = x * scale_factor;
y_scaled = y * scale_factor;

% Plot the Scaled 3D Track
figure;
plot3(x_scaled, y_scaled, z, 'k-', 'LineWidth', 2); % Plot the 3D track
grid on;
title('3D Track Visualization (Scaled)');
xlabel(['X (meters) scaled by ', num2str(scale_factor)]);
ylabel(['Y (meters) scaled by ', num2str(scale_factor)]);
zlabel('Elevation (meters)');

set(gca, 'XDir', 'reverse'); % Reverse the direction of the X-axis

% Set Axis Limits
xlim([min(x_scaled), max(x_scaled)]);
ylim([min(y_scaled), max(y_scaled)]);
zlim([min(z), max(z)]);
axis equal; % Maintain aspect ratio
rotate3d on; % Enable interactive rotation



