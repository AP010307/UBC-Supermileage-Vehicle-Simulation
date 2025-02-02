% Load track coordinates from a CSV file
track_data = readmatrix('sem_2023_us.csv'); % Replace with your file name

% Extract X, Y, and Z (if applicable)
x = track_data(:, 1); % First column: X-coordinates
y = track_data(:, 2); % Second column: Y-coordinates
z = track_data(:, 3); % Third column: Z-coordinates

% 2D Track Visualization
figure;
plot(x, y, 'k-', 'LineWidth', 2); % 2D track plot
grid on;
title('2D Track Visualization');
xlabel('X (meters)');
ylabel('Y (meters)');
axis equal; % Maintain aspect ratio

% 3D Track Visualization
figure;
plot3(x, y, z, 'k-', 'LineWidth', 2); % 3D track plot
grid on;
title('3D Track Visualization');
xlabel('X (meters)');
ylabel('Y (meters)');
zlabel('Elevation (meters)');
axis equal; % Maintain aspect ratio

%2D view
hold on;
plot(x(1), y(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); % Start point
plot(x(end), y(end), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % End point
legend('Track', 'Start', 'Finish');

%3D view
hold on;
plot3(x(1), y(1), z(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); % Start point
plot3(x(end), y(end), z(end), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % End point
legend('Track', 'Start', 'Finish');


figure;
scatter3(x, y, z, 15, z, 'filled'); % Points colored by elevation
colorbar; % Add color scale
title('3D Track with Elevation Colors');
xlabel('X (meters)');
ylabel('Y (meters)');
zlabel('Elevation (meters)');
axis equal;

%Animate
figure;
plot3(x, y, z, 'k-', 'LineWidth', 2); % Plot the track
hold on;
car_plot = plot3(x(1), y(1), z(1), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % Initial position

% Animation loop
for i = 1:length(x)
    set(car_plot, 'XData', x(i), 'YData', y(i), 'ZData', z(i)); % Update position
    pause(0.05); % Adjust pause time for animation speed
end

% 2D Visualization
figure;
subplot(1, 2, 1);
plot(x, y, 'k-', 'LineWidth', 2);
hold on;
plot(x(1), y(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); % Start
plot(x(end), y(end), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Finish
title('2D Track Visualization');
xlabel('X (meters)');
ylabel('Y (meters)');
axis equal;
legend('Track', 'Start', 'Finish');

% 3D Visualization
subplot(1, 2, 2);
plot3(x, y, z, 'k-', 'LineWidth', 2);
hold on;
plot3(x(1), y(1), z(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); % Start
plot3(x(end), y(end), z(end), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Finish
grid on;
title('3D Track Visualization');
xlabel('X (meters)');
ylabel('Y (meters)');
zlabel('Elevation (meters)');
axis equal;
legend('Track', 'Start', 'Finish');
