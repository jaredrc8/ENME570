% ENME 570 Lab

clc
clear
close all

% Load the data
data = load('naca0012.dat');
cl = load('Cl.txt');

% Top of the airfoil
x = data(1:length(data)/2+1, 1);
y = data(1:length(data)/2+1, 2);

% Fit a spline to the data
pp = spline(x, y);

% Plot the original data points
figure;
hold on;
grid on;
plot(x, y, 'o', 'DisplayName', 'Original Data');
plot(x, -y,'o', 'DisplayName', 'Original Data');
axis([0 1 -1 1])

% Plot the fitted spline
fitted_x = linspace(min(x), max(x), 1000);
fitted_y = ppval(pp, fitted_x);
plot(fitted_x, fitted_y, '-', 'DisplayName', 'Fitted Spline');
plot(fitted_x, -fitted_y, '-', 'DisplayName', 'Fitted Spline');

% Usage:
% x_query = 0.8;  % Distance from leading edge along the chord
% normal_vec = normalAtX(x_query, pp);

% Plot the normal vector at x_query
% normal_start = [x_query, ppval(pp, x_query)];
% quiver(normal_start(1), normal_start(2), normal_vec(1), normal_vec(2), 0.1, 'r', 'LineWidth', 1, 'DisplayName', 'Normal Vector');

% Display the results
% fprintf('The normal vector at x = %.2f is: [%.4f, %.4f]\n', x_query, normal_vec(1), normal_vec(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify chord length and pressure tap locations (mm)
chord = 152;
tap_locations = [0.76; 1.52; 3.81; 7.62; 11.43; 15.24; 19.05; 22.86; 38; 41.15; 62; 59.44; 80.77; 77.73; 101.35; 96.02; 121.92; 114.30; 137.16; 129.54];

% For every pressure tap, calculate the normal vector
for i = 1:length(tap_locations)
    normal_vec = normalAtX(tap_locations(i,1)/chord , pp);
    tap_vectors(i,1) = normal_vec(1);
    % Every odd numbered tap is for the lower airfoil surface (-y vector)
    if mod(i,2) == 0
        tap_vectors(i,2) = - normal_vec(2);
        tap_locations(i,2) = ppval(pp, tap_locations(i,1)/chord);
        quiver(tap_locations(i,1)/chord, -ppval(pp, tap_locations(i,1)/chord), tap_vectors(i,1)/5, tap_vectors(i,2)/5);
    else
        tap_vectors(i,2) = normal_vec(2);
        tap_locations(i,2) = - ppval(pp, tap_locations(i,1)/chord);
        quiver(tap_locations(i,1)/chord, ppval(pp, tap_locations(i,1)/chord), tap_vectors(i,1)/5, tap_vectors(i,2)/5);
    end
end

panel_lengths(1,1) = pdist2([tap_locations(1,1), tap_locations(1,2)] , [tap_locations(2,1), tap_locations(2,2)] , 'euclidean');
for i = 2:length(tap_locations(:,1))-1
    panel_lengths(i,1) = pdist2([tap_locations(i,1), tap_locations(i,2)] , [tap_locations(i+1,1), tap_locations(i+1,2)] , 'euclidean')/2;
    
end

% From Cl import, sum all force vectors for each angle of attack
for i = 1:length(cl(1,:)) % for each angle of attack (by column in cl)
    lift(i) = 0;
    drag(i) = 0;
    for j = 1:length(cl(:,1)) %for each pressure tap per angle of attack
        lift(i) = lift(i) + cl(j,i) * tap_vectors(j,2);
        drag(i) = drag(i) + cl(j,i) * tap_vectors(j,1);
    end
end

% Define a function to calculate the normal vector at a given x position
function normal_vector = normalAtX(x_val, pp)
    % Calculate the first derivative (slope)
    dydx = ppval(fnder(pp, 1), x_val);
    
    % Tangent vector (normalized)
    tangent_vector = [1, dydx] / norm([1, dydx]);
    
    % Normal vector (90 degrees rotated)
    normal_vector = [-tangent_vector(2), tangent_vector(1)];
end