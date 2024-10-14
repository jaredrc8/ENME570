% NACA0012 Aerofoil Analysis
clear
clc
close all

% Read Data
Data_force = readmatrix('NACA0012_Lift.xlsx', 'range', 'b6:b1103');
Data_Angle = readmatrix('NACA0012_Lift.xlsx', 'range', 'au6:au1103');

% Constants
rho_water = 997; % Density of water [kg/m^3]
rho_oil_specific = 0.87; % Specific density of pitot tube manometer
rho_oil = rho_water * rho_oil_specific; % Density of manometer fluid [kg/m^3]
rho_air = 1.06; % Density of air [kg/m^3]
g = 9.81; % Gravity [m/s^2]

% Pitot Tube Reading
p_tube = 40; % Pitot tube reading [mm]
U_inf = sqrt((2 * rho_oil * g * p_tube / 1000) / rho_air); % Free Stream Velocity [m/s]

% Wing Dimensions
s = 0.300;   % Wing span [m]
c = 0.15001; % Chord [m]
A = s * c;   % Reference Area [m^2]

% Unique angles and their corresponding mean lift forces
alpha = unique(Data_Angle); % Get unique angles of attack
mean_lift_forces = zeros(size(alpha)); % Preallocate array for mean lift forces

% Calculate Mean Lift Force for each Angle of Attack
for i = 1:length(alpha)
    % Find indices of measurements corresponding to the current angle
    indices = (Data_Angle == alpha(i));
    
    % Calculate mean lift force for these indices
    mean_lift_forces(i) = mean(Data_force(indices));
end

mean_lift_forces = abs( mean_lift_forces ) - mean_lift_forces(1);

% Calculate Lift Coefficient
C_L = mean_lift_forces / (0.5 * rho_air * U_inf^2 * A); % Lift Coefficient

% Plotting CL Force vs Angle of Attack
figure;
plot(alpha, C_L, 'b-o'); % Plot CL vs alpha
xlabel('Angle of Attack (degrees)');
ylabel('CL');
title('CL vs Angle of Attack for NACA0012 Airfoil');
grid on;


% Read Data For Drag
Data_force_2 = readmatrix('NACA0012_Drag.xlsx', 'range', 'b6:b1102');
Data_Angle_2 = readmatrix('NACA0012_Drag.xlsx', 'range', 'au6:au1102');

% Unique angles and their corresponding mean lift forces
alpha_1 = unique(Data_Angle_2); % Get unique angles of attack
mean_drag_forces = zeros(size(alpha_1)); % Preallocate array for mean lift forces

% Calculate Mean Lift Force for each Angle of Attack
for i = 1:length(alpha_1)
    % Find indices of measurements corresponding to the current angle
    indices_1 = (Data_Angle_2 == alpha_1(i));
    
    % Calculate mean lift force for these indices
    mean_drag_forces(i) = mean(Data_force_2(indices_1));
end


% Calculate drag Coefficient
C_D = mean_drag_forces / (0.5 * rho_air * U_inf^2 * A); % Lift Coefficient

% Plotting Cd vs Angle of Attack
figure;
plot(alpha_1, C_D, 'b-o'); % Plot Cd vs alpha
xlabel('Angle of Attack (degrees)');
ylabel('CD');
title('CD vs Angle of Attack for NACA0012 Airfoil');
grid on;

% Plotting Cd vs CL
figure;
plot(C_D,C_L,'b-o'); % Plot Cd vs CL
xlabel('CD');
ylabel('CL');
title('CL vs CD for NACA0012 Airfoil');
grid on;