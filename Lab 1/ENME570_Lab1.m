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
q = 0.5 * air_density * airspeed^2;

% Load the data
data = load('naca0012.dat');
pressure = load('pressure.txt');

% Top of the airfoil
x = data(1:length(data)/2+1, 1);
y = data(1:length(data)/2+1, 2);

% Fit a spline to the data
pp = spline(x, y);

% Plot the original data points
figure;
hold on;
grid on;
% plot(x, y, 'o', 'DisplayName', 'Original Data');
% plot(x, -y,'o', 'DisplayName', 'Original Data');
axis([0 1 -1 1])

% Plot the fitted spline
fitted_x = linspace(min(x), max(x), 1000);
fitted_y = ppval(pp, fitted_x);
plot(fitted_x, fitted_y, '-', 'DisplayName', 'Fitted Spline');
plot(fitted_x, -fitted_y, '-', 'DisplayName', 'Fitted Spline');

%%

% For every pressure tap, calculate the normal vector
for i = 1:length(tap_locations)
    normal_vec = normalAtX(tap_locations(i,1)/chord , pp);
    tap_vectors(i,1) = normal_vec(1);
    % Every even numbered tap is for the lower airfoil surface (-y vector)
    if mod(i,2) == 0
        tap_vectors(i,2) = - normal_vec(2);
        tap_locations(i,2) = - chord * ppval(pp, tap_locations(i,1)/chord);
        quiver(tap_locations(i,1)/chord, -ppval(pp, tap_locations(i,1)/chord), tap_vectors(i,1)/5, tap_vectors(i,2)/5);
    else
        tap_vectors(i,2) = normal_vec(2);
        tap_locations(i,2) = chord * ppval(pp, tap_locations(i,1)/chord);
        quiver(tap_locations(i,1)/chord, ppval(pp, tap_locations(i,1)/chord), tap_vectors(i,1)/5, tap_vectors(i,2)/5);
    end
end

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

% Lift Calculation
% For each angle of attack
for AoA = 1:length(pressure(1,:))
    % For each panel/pressure tap
    for i = 1:length(tap_locations(:,1))
        dA = (span * panel_lengths(i,1)) / 10^6;
        % Direction of tap_vector accounts for upper/lower surface
        dL(i,AoA) = pressure(i,AoA) * dA * tap_vectors(i,2);
    end
    lift(AoA) = sum(dL(:,AoA)) - sum(dL(:,1));
end

% Coefficient of Lift Calculation
for i = 1:length(lift)
    CL(i) = lift(i)/(q * (span * chord) / 10^6);
end

figure
hold on
title("CL v AoA of NACA0012 Airfoil (Pressure Distribution)")
plot([0,2,4,6,8,10,12,14,16],CL)


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