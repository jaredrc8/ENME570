% ENME 570 Lab

clc
clear
close all

% Specify properties and pressure tap locations (mm)
chord = 152;
tap_locations = [0.76; 1.52; 3.81; 7.62; 11.43; 15.24; 19.05; 22.86; 38; 41.15; 62; 59.44; 80.77; 77.73; 101.35; 96.02; 121.92; 114.30; 137.16; 129.54];
span = 300;
air_density = 1.06; % kg/m3
air_nu = 0.00001506; % m2/s
water_density = 997; % kg/m3
oil_density = 867.4; % kg/m3
g = 9.81; % m/s2
airspeed = sqrt((2*oil_density*g*40/1000)/(air_density));
q = 0.5 * air_density * airspeed^2; % dynamic pressure


% Load the data
data = load('naca0012.dat');
pressure = readmatrix("AirfoilLab_Template.xlsx","range","l19:t39");
pressure_literature = readmatrix("wpd_datasets.csv");

% Top of the airfoil
x = data(1:length(data)/2+1, 1);
y = data(1:length(data)/2+1, 2);

% Fit a spline to the data
pp = spline(x, y);

% Plot the original data points
figure;
hold on;
grid on;
title("NACA0012 Airfoil Visual Representation with Pressure Taps")
axis([-0.2 1.1 -0.5 0.5])
% plot(x, y, 'or', 'DisplayName', 'Original Data');
% plot(x, -y,'or', 'DisplayName', 'Original Data');

% Plot the fitted spline
fitted_x = linspace(min(x), max(x), 1000);
fitted_y = ppval(pp, fitted_x);
plot(fitted_x, fitted_y, 'r', 'DisplayName', 'Airfoil');
plot(fitted_x, -fitted_y, 'r', 'DisplayName', 'Airfoil');

%%

% For every pressure tap, calculate the normal vector
for i = 1:length(tap_locations)
    normal_vec = normalAtX(tap_locations(i,1)/chord , pp);
    tap_vectors(i,1) = normal_vec(1);
    % Every even numbered tap is for the lower airfoil surface (-y vector)
    if mod(i,2) == 0
        tap_vectors(i,2) = - normal_vec(2);
        tap_locations(i,2) = - chord * ppval(pp, tap_locations(i,1)/chord);
        quiver(tap_locations(i,1)/chord, -ppval(pp, tap_locations(i,1)/chord), tap_vectors(i,1), tap_vectors(i,2),0.2,'b','DisplayName','Pressure Taps');
    else
        tap_vectors(i,2) = normal_vec(2);
        tap_locations(i,2) = chord * ppval(pp, tap_locations(i,1)/chord);
        quiver(tap_locations(i,1)/chord, ppval(pp, tap_locations(i,1)/chord), tap_vectors(i,1), tap_vectors(i,2),0.2,'b','DisplayName','Pressure Taps');
    end
end
legend('Airfoil', '','Pressure Taps');

% UNFINISHED SECTION: Visual representation of pressure tap locations and
% panels for each pressure tap.
% 
% 
% panel_points = nan(length(tap_locations(:,1)),2);
% panel_points(1,1) = 0;
% panel_points(1,2) = 0;
% panel_points(length(tap_locations(:,1)),1) = chord;
% panel_points(length(tap_locations(:,1)),2) = 0;
% for i = length(tap_locations(:,1))-1:-1:3
%     panel_points(i,1) = 0.5 * (tap_locations(i,1) + tap_locations(i-2,1));
%     panel_points(i,2) = 0.5 * (tap_locations(i,2) + tap_locations(i-2,2));
% end
% panel_points(2,1) = 0.5 * (tap_locations(i,1) + tap_locations(i-1,1));
% panel_points(2,2) = 0.5 * (tap_locations(i,2) + tap_locations(i-1,2));
% 
% figure
% hold on
% grid on
% plot(chord*fitted_x, chord*fitted_y, '-', 'DisplayName', 'Fitted Spline');
% plot(chord*fitted_x, -chord*fitted_y, '-', 'DisplayName', 'Fitted Spline');
% plot(panel_points(:,1),panel_points(:,2),'+','DisplayName','Panel Points');
% plot(tap_locations(:,1),tap_locations(:,2),'o','DisplayName','Tap Locations');
% hold off


% Calculate lengths of panels associated with each pressure tap
% Assume even number of pressure taps
% First panels extend to leading edge (offsets pressure tap from midpoint, index 1,2)
panel_lengths(1,1) = pdist2(chord*[0,0],[tap_locations(1,1),tap_locations(1,2)],"euclidean") + 0.5*pdist2([tap_locations(1,1), tap_locations(1,2)] , [tap_locations(3,1), tap_locations(3,2)] , 'euclidean');
panel_lengths(2,1) = pdist2(chord*[0,0],[tap_locations(2,1),tap_locations(2,2)],"euclidean") + 0.5*pdist2([tap_locations(2,1), tap_locations(2,2)] , [tap_locations(4,1), tap_locations(4,2)] , 'euclidean');
% Last panels extends to trailing edge (offsets pressure tap from midpoint, index 19,20)
panel_lengths(length(tap_locations(:,1)),1) = pdist2(chord*[1,0],[tap_locations(end,1),tap_locations(end,2)],"euclidean") + 0.5*pdist2([tap_locations(end,1), tap_locations(end,2)] , [tap_locations(end-2,1), tap_locations(end-2,2)] , 'euclidean');
panel_lengths(length(tap_locations(:,1))-1,1) = pdist2(chord*[1,0],[tap_locations(end-1,1),tap_locations(end-1,2)],"euclidean") + 0.5*pdist2([tap_locations(end-1,1), tap_locations(end-1,2)] , [tap_locations(end-3,1), tap_locations(end-3,2)],"euclidean");
% Subsequent panels contain pressure tap at exact midpoint
% Top of airfoil (odd numbers)
for i = length(tap_locations(:,1))-2:-2:4
    panel_lengths(i,1) = 0.5*pdist2([tap_locations(i,1),tap_locations(i,2)],[tap_locations(i+2,1),tap_locations(i+2,2)],'euclidean') + 0.5*pdist2([tap_locations(i,1),tap_locations(i,2)],[tap_locations(i-2,1),tap_locations(i-2,2)],'euclidean');
end
% Bottom of airfoil (even numbers)
for i = length(tap_locations(:,1))-3:-2:3
    panel_lengths(i,1) = 0.5*pdist2([tap_locations(i,1),tap_locations(i,2)],[tap_locations(i+1,1),tap_locations(i+1,2)],"euclidean") + 0.5*pdist2([tap_locations(i-1,1),tap_locations(i-1,2)],[tap_locations(i,1),tap_locations(i,2)],"euclidean");
end
% Check if sum of panel lengths approximates total airfoil circumference
airfoil_circumference = 2 * chord * integral(@(t) sqrt(1 + ppval(fnder(pp), t).^2), min(x), max(x));
disp("Sum of Panel Lengths: " + sum(panel_lengths) + "mm");
disp("Airfoil Circumference: " + airfoil_circumference + "mm");
disp("Discretization Error: " + ((sum(panel_lengths)-airfoil_circumference)/airfoil_circumference*100) + "%");

% Lift/Drag Calculation
static_pressure = pressure(end,:)';
% For each angle of attack
for AoA = 1:length(pressure(1,:))
    % For each panel/pressure tap
    for i = 1:length(tap_locations(:,1))
        dA = (span * panel_lengths(i,1)) / 10^6;
        % Direction of tap_vector accounts for upper/lower surface
        dL(i,AoA) = (pressure(i,AoA) - static_pressure(AoA,1)) * dA * tap_vectors(i,2);
        dD(i,AoA) = (pressure(i,AoA) - static_pressure(AoA,1)) * dA * tap_vectors(i,1);
        Cp(i,AoA) = -(pressure(i,AoA) - static_pressure(AoA,1)) / q;
    end
    % Sum differential lift to equal total lift, subtracted by lift at
    % AoA=0 as a correction
    lift(AoA) = (sum(dL(:,AoA)) - sum(dL(:,1))) * cos(deg2rad(AoA*2-2));
    drag(AoA) = (sum(dD(:,AoA)) - sum(dD(:,1))) * cos(deg2rad(AoA*2-2));
end

% Coefficient Calculations
for i = 1:length(lift)
    CL(i) = lift(i)/(q * (span * chord) / 10^6);
    CD(i) = drag(i)/(q * (span * chord) / 10^6);
end

figure
hold on
grid on
title("Cl v AoA of NACA0012 Airfoil (Pressure Distribution)")
plot(linspace(0,16,9),CL,'b-o')
xlabel("Angle of Attack (degrees)")
ylabel("Cl")

figure
hold on
grid on
title("Cp Distribution of NACA0012")
Cp_upper = Cp(1:2:19,:);
Cp_lower = Cp(2:2:20,:);
% AoA = 2
plot(tap_locations(2:2:20,1)/chord,Cp_upper(:,2),'r-o','DisplayName','Upper Airfoil Surface');
plot(tap_locations(1:2:19,1)/chord,Cp_lower(:,2),'r-+', 'DisplayName', 'Lower Airfoil Surface');
% AoA = 6
plot(tap_locations(2:2:20,1)/chord,Cp_upper(:,4),'b-o','DisplayName','Upper Airfoil Surface');
plot(tap_locations(1:2:19,1)/chord,Cp_lower(:,4),'b-+', 'DisplayName', 'Lower Airfoil Surface');
% AoA = 10
plot(tap_locations(2:2:20,1)/chord,Cp_upper(:,6),'k-o','DisplayName','Upper Airfoil Surface');
plot(tap_locations(1:2:19,1)/chord,Cp_lower(:,6),'k-+', 'DisplayName', 'Lower Airfoil Surface');
% AoA = 14
plot(tap_locations(2:2:20,1)/chord,Cp_upper(:,8),'m-o','DisplayName','Upper Airfoil Surface');
plot(tap_locations(1:2:19,1)/chord,Cp_lower(:,8),'m-+', 'DisplayName', 'Lower Airfoil Surface');
legend('Upper Airfoil Surface a=2','Lower Airfoil Surface a=2','Upper Airfoil Surface a=6','Lower Airfoil Surface a=6','Upper Airfoil Surface a=10','Lower Airfoil Surface a=10','Upper Airfoil Surface a=14','Lower Airfoil Surface a=14','Location','best');
xlabel("Nondimensional position, x/chord")
ylabel("Cp")
axis([0 0.92 -5 2])

figure
hold on
grid on
title("Cp Distribution of NACA0012 at AoA=10");
plot(pressure_literature(:,1),pressure_literature(:,2),'r-','DisplayName','Literature');
plot(pressure_literature(:,3),pressure_literature(:,4),'r-','DisplayName','Literature');
plot(tap_locations(2:2:20,1)/chord,Cp_upper(:,6),'b-o','DisplayName','Experimental');
plot(tap_locations(1:2:19,1)/chord,Cp_lower(:,6),'b-o','DisplayName','Experimental');
axis([0 0.92 -5 2])
xlabel("Nondimensional position, x/chord")
ylabel("Cp")
legend('Literature','','','Experimental','Location','best')

%%
% FUNCTIONS

% Define a function to calculate the normal vector at a given x position
function normal_vector = normalAtX(x_val, pp)
    % Calculate the first derivative (slope)
    dydx = ppval(fnder(pp, 1), x_val);

    % Tangent vector (normalized)
    tangent_vector = [1, dydx] / norm([1, dydx]);
    
    % Normal vector (90 degrees rotated)
    normal_vector = [-tangent_vector(2), tangent_vector(1)];
end