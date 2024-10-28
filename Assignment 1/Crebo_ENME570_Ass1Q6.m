% Jared Crebo 30085839
% ENME 570 Assignment 1 Question 6

clc
clear
close all

%% Question 6a
% Adjust figure settings
axis_size = 2;

% i)
% Flow Dominated by Freestream Velocity
velocity_field(10, 2, axis_size); % U = 10m/s, Lambda = 2m^2/s
title("6a i) Velocity Field of Flow Dominated by Freestream Velocity")
hold off
% Flow Dominated by Vortex Strength
velocity_field(2, 10, axis_size); % U = 2m/s, Lambda = 10m^2/s
title("6a i) Velocity Field of Flow Dominated by Vortex Strength")
hold off

% ii)
% Velocity Field with Circle Contour (d = 0.6m)
velocity_field(5, 5, axis_size); % U = 5m/s, Lambda = 5m^2/s
title("6a ii) Velocity Field with Circle Contour (d = 0.6m)");
circulation("circle", 0.6, 0, 0); % circle contour, 0.6m diameter
hold off

% iii)
% Velocity Field with Circle Contour (d = 1m)
velocity_field(5, 5, axis_size); % U = 5m/s, Lambda = 5m^2/s
title("6a iii) Velocity Field with Circle Contour (d = 1m)")
circulation("circle", 1, 0, 0); % circle contour, 1m diameter
hold off
% Velocity Field with Square Contour (l = 1m)
velocity_field(5, 5, axis_size); % U = 5m/s, Lambda = 5m^2/s
title("6a iii) Velocity Field with Square Contour (l = 1m)")
circulation("square", 1, 0, 0); % square contour, 1m length
hold off
% Velocity Field with Ellipse Contour (a = 1m)
velocity_field(5, 5, axis_size); % U = 5m/s, Lambda = 5m^2/s
title("6a iii) Velocity Field with Ellipse Contour (a = 1m)")
circulation("ellipse", 1, 0, 0); % ellipse contour, 1m semimajor axis
hold off

% iv)
% Velocity Field with Offset Circle Contour
velocity_field(5, 5, axis_size);
title("6a iv) Velocity Field with Offset Circle Contour");
circulation("circle", 0.6, 1, 0);
circulation("circle", 1, -1, -0.5);
hold off

%% Question 6b
% Adjust figure settings
axis_size = 3;

% ii), iii)
% Velocity Field Around Two Vortices
velocity_field2(5, 5, axis_size, 2);
title("6b ii), iii) Velocity Field Around Two Vortices");
circulation("ellipse", 2, 1, 0);
hold off

%% Question 6c
% Adjust figure settings
axis_size = 5;

% i)
% Velocity Field Around Three Vortices
velocity_field3(5, 5, axis_size, 1.5);
title("6c i) Velocity Field Around Three Vortices");
circulation("ellipse", 3, 1.5, 0)
hold off

%% Question 6d
% Adjust figure settings
axis_size = 10;

% i)
% Velocity Field Around N Vortices
velocity_fieldn(5, 5, axis_size, 2, 5); % U=5, Lambda=50, axis=10, a=2, n=5
title("6d i) Velocity Field Around N=5 Vortices");
circulation("ellipse", 5.5, 4, 0); % semimajor axis = 5.5, x = 4, y = 0
circulation("square", 10, 4, 0);
circulation("circle", 10, 4, 0);
hold off

%% Question 6e (Bonus)
% Define thin airfoil
chord = linspace(0, 1, 100); % LE to TE normalized, discretized
[x, y] = meshgrid(linspace(-0.1,1.1,100), linspace(-0.1,0.1,100));
alpha = deg2rad(linspace(0, 12, 4));
U = 5; % freestream velocity m/s

% For each angle of attack
for i = 1:length(alpha)
    theta = acos(1 - 2*chord); % theta of chord discretization
    gamma = 2 * alpha(i) * U * ((1 + cos(theta))/sin(theta)); % vortex strength
    r_sq = (x - chord(i)).^2 + y.^2; % r^2
    
    u = U + zeros(size(x));
    v = zeros(size(y));

    for j = 1:length(chord)-1
        % Compute induced velocity at each point along chord
        dx = diff(x);
        dy = diff(y);
        u_induced = - gamma * y ./ (2 * pi .* r_sq) * dx(j) ;
        v_induced = gamma * (x - chord(i)) ./ (2 * pi .* r_sq) * dy(j) ;
        u = u + u_induced;
        v = v + v_induced;
    end

    figure;
    hold on
    velocity = quiver(x,y,u,v,'Color','m','DisplayName','Velocity Vector Field');
    airfoil = plot(chord, linspace(0,0,length(chord)), 'k', 'LineWidth', 2,'DisplayName','Thin Airfoil');
    streamlines = streamline(x,y,u,v,x(1,1), y(1,1));
    set(streamlines, 'DisplayName', 'Streamline');
    for k = 3:2:length(y(:,1))
        streamline(x,y,u,v,x(1,1), y(k,1));
    end
    xlabel('x');
    ylabel('y');
    title(sprintf('Velocity Field around Thin Airfoil at $\\alpha = %1.f^\\circ$', rad2deg(alpha(i))), 'Interpreter', 'latex');
    axis([-0.1 1.1 -0.1 0.1]);
    legend([velocity,airfoil, streamlines]);
    hold off
end

%% Functions

% Create a velocity field(freestream velocity (m/s), vortex strength (m2/s)
function velocity_field(U, Lambda, axis_size)
    % Create a grid
    global x y u v
    
    [x, y] = meshgrid(-axis_size:axis_size/40:axis_size, -axis_size:axis_size/40:axis_size);
    
    % Calculate velocity components
    u = U - Lambda/(2*pi) .* y ./ (x.^2 + y.^2);
    v = Lambda/(2*pi) .* x ./ (x.^2 + y.^2);

    % Plot the velocity field
    figure;
    hold on
    quiver(x, y, u, v);
    axis([-axis_size axis_size -axis_size axis_size]);
    axis square;
    title('Velocity Field');
    xlabel('x');
    ylabel('y');
    for i = -axis_size:axis_size/30:axis_size
        set(streamline(x, y, u, v, -axis_size, i), 'DisplayName', 'Streamline');
    end
    legend('Velocity Vector','Streamline');
end

% Create a velocity field and define 2 vortices separated by distance "a"
function velocity_field2(U, Lambda, axis_size, a)
    global x y u v
    
    [x, y] = meshgrid(-axis_size:axis_size/40:axis_size, -axis_size:axis_size/40:axis_size);
    
    % Calculate velocity components
    u = U - Lambda/(4*pi) .* y .* (1./(x.^2 + y.^2) + 1./((x - a).^2 + y.^2));
    v = Lambda/(4*pi) .* (x./(x.^2 + y.^2) + (x - a)./((x - a).^2 + y.^2));

    % Plot the velocity field
    figure;
    hold on
    quiver(x, y, u, v);
    axis([-axis_size axis_size -axis_size axis_size]);
    axis square;
    title('Velocity Field');
    xlabel('x');
    ylabel('y');
    for i = -axis_size:axis_size/30:axis_size
        set(streamline(x, y, u, v, -axis_size, i), 'DisplayName', 'Streamline');
    end
    legend('Velocity Vector','Streamline');
end

% Create a velocity field and define 3 vortices separated by distance "a"
function velocity_field3(U, Lambda, axis_size, a)
    global x y u v
    
    [x, y] = meshgrid(-axis_size:axis_size/40:axis_size, -axis_size:axis_size/40:axis_size);
    
    % Calculate velocity components
    u = U - Lambda/(6*pi) .* y .* (1./(x.^2 + y.^2) + 1./((x - a).^2 + y.^2) + 1./((x - 2*a).^2 + y.^2));
    v = Lambda/(6*pi) .* (x./(x.^2 + y.^2) + (x - a)./((x - a).^2 + y.^2) + (x - 2*a)./((x - 2*a).^2 + y.^2));

    % Plot the velocity field
    figure;
    hold on
    quiver(x, y, u, v);
    axis([-axis_size axis_size -axis_size axis_size]);
    axis square;
    title('Velocity Field');
    xlabel('x');
    ylabel('y');
    for i = -axis_size:axis_size/20:axis_size
        set(streamline(x, y, u, v, -axis_size, i), 'DisplayName', 'Streamline');
    end
    legend('Velocity Vector','Streamline');
end

% Create a velocity field and define n vortices separated by distance "a"
function velocity_fieldn(U, Lambda, axis_size, a, n)
    global x y u v
    
    [x, y] = meshgrid(-axis_size:axis_size/40:axis_size, -axis_size:axis_size/40:axis_size);
    n = 1:1:n;

    % Calculate velocity components
    u_ = 0;
    v_ = 0;
    for i = n
        u_ = u_ + 1./((x - (n(i)-1)*a).^2 + y.^2);
        v_ = v_ + (x - (n(i)-1)*a) ./ ((x - (n(i)-1)*a).^2 + y.^2);
    end

    u = U - Lambda.*y./(n(end)*2*pi) .* u_;
    v = Lambda./(n(end)*2*pi) .* v_;
  
    % Plot the velocity field
    figure;
    hold on
    quiver(x, y, u, v);
    axis([-axis_size axis_size -axis_size axis_size]);
    axis square;
    title('Velocity Field');
    xlabel('x');
    ylabel('y');
    for i = -axis_size:axis_size/20:axis_size
        set(streamline(x, y, u, v, -axis_size, i), 'DisplayName', 'Streamline');
    end
    legend('Velocity Vector','Streamline');
end


% Shape: square, circle, ellipse
function circulation = circulation(shape, diameter, x_o, y_o)
    global x y u v
    if shape == "square"
        % Define the square
        L = diameter;
        
        x1 = x_o + linspace(-L/2, L/2, 100);  % Bottom edge (y = -L/2)
        y1 = y_o + -L/2 * ones(1, 100);
        
        x2 = x_o + L/2 * ones(1, 100);        % Right edge (x = L/2)
        y2 = y_o + linspace(-L/2, L/2, 100);
        
        x3 = x_o + linspace(L/2, -L/2, 100);  % Top edge (y = L/2)
        y3 = y_o + L/2 * ones(1, 100);
        
        x4 = x_o + -L/2 * ones(1, 100);       % Left edge (x = -L/2)
        y4 = y_o + linspace(L/2, -L/2, 100);
        
        % Concatenate the points to form the full contour
        x_square = [x1, x2, x3, x4];
        y_square = [y1, y2, y3, y4];
        
        % Plot square shape
        plot(x_square,y_square,'r','LineWidth', 1, 'DisplayName', 'Circulation Contour');

        dx = diff(x_square);
        dy = diff(y_square);

        u_square = interp2(x, y, u, x_square, y_square);
        v_square = interp2(x, y, v, x_square, y_square);
        
        % Calculate circulation about the contour, line integral summation
        circulation = sum(u_square(1:end-1) .* dx + v_square(1:end-1) .* dy);
    elseif shape == "circle"
        % Define a circular path around the vortex
        theta = linspace(0, 2*pi, 100);
        radius = diameter / 2;
        cx = x_o + radius * cos(theta);
        cy = y_o + radius * sin(theta);

        % Interpolate velocities at the circle points
        u_circle = interp2(x, y, u, cx, cy);
        v_circle = interp2(x, y, v, cx, cy);

        % Plot circulation shape
        plot(cx,cy,'r','LineWidth', 1, 'DisplayName', 'Circulation Contour')
        
        % Calculate circulation about the contour, line integral summation
        ds = 2 * pi * radius / length(theta); % arc length element
        vt = u_circle .* (-sin(theta)) + v_circle .* cos(theta); % tangential velocity around contour
        circulation = sum(vt * ds);
    elseif shape == "ellipse"
        % Define the ellipse
        a = diameter; % semi major axis
        b = diameter / 2; % semi minor axis
        theta = linspace(0, 2*pi, 100);
        x_ellipse = x_o + a * cos(theta);
        y_ellipse = y_o + b * sin(theta);

        % Interpolate velocities at ellipse
        u_ellipse = interp2(x, y, u, x_ellipse, y_ellipse);
        v_ellipse = interp2(x, y, v, x_ellipse, y_ellipse);

        % Plot ellipse
        plot(x_ellipse,y_ellipse,'r','LineWidth', 1, 'DisplayName', 'Circulation Contour');

        dx = diff(x_ellipse);
        dy = diff(y_ellipse);
        
        % Calculate circulation about the contour, line integral summation
        circulation = sum(u_ellipse(1:end-1) .* dx + v_ellipse(1:end-1) .* dy);
    else
        error("Not configured for this shape type. Check your spelling.");
    end
    disp("Circulation: " + num2str(circulation));
end