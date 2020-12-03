% Extent-Compatible Control Barrier Functions
% NOTE: Please run the "init.m" file from the Robotarium package prior to
%       executing this code
% Authors: Mohit Srinivasan, Matt Abate, Gustav Nilsson, and Samuel Coogan
% Date Modified: 03/25/2020

clear all;
close all;
clc;

% Parameters 
% The domain
x1_dom = -1.6:0.01:1.6; x2_dom = -1.6:0.01:1.6;
y1_dom = -1.6:0.01:1.6; y2_dom = -1.6:0.01:1.6;

% The safe set the robot should stay within
a = 1;
b = 0.8;
P_safe = [1/a^2 0; 0 1/b^2];

% The extent set for the robot
%shape = 0.1*[1 0; 0 1.333];
shape = 0.1*[1 0; 0 2];

% Number of robots
N = 1;

% Number of iterations 
iterations = 1000;

% Compute the reduced safeset
extent_max_radius = max(max(shape));
a_r = a - extent_max_radius;
b_r = b - extent_max_radius;
P_safe_r = [1/a_r^2 0; 0 1/b_r^2];

% Which controller we want to use
cont = 'extent'; % Pick 'extent', 'sos', or 'point'

%% Setup for the SAMPLE based controller
% Number of samples
if strcmp(cont,'extent')
    num_samp = 2000;

    % Gamma paramter for the sampling based
    gamma = 0.06;

    % Determine the sampling distance (tau)
    th = 0:(2*pi/num_samp):2*pi;
    latest_sampling_point = [sqrt(abs(cos(th(1))))*shape(1, 1)*sign(cos(th(1))); ...
        sqrt(abs(sin(th(1))))*shape(2, 2)*sign(sin(th(1)))];
    tau = inf;
    for i=2:length(th)
        new_sampling_point = [sqrt(abs(cos(th(i))))*shape(1, 1)*sign(cos(th(i))); ...
            sqrt(abs(sin(th(i))))*shape(2, 2)*sign(sin(th(i)))];

        if norm(latest_sampling_point-new_sampling_point) < tau
            tau = norm(latest_sampling_point-new_sampling_point);
        end
        latest_sampling_point = new_sampling_point;
    end

    % Determine the constants
    X_dom = [x1_dom; x2_dom];
    Y_dom = [y1_dom; y2_dom];
    dh_dx = 2*X_dom'*P_safe;
    dE2_dxdy = [];
    for i = 1:length(x1_dom)
        for j = 1:length(y1_dom)
            temp = [2*norm([2*(x1_dom(i)-y1_dom(j)) 0; 0 2*(x2_dom(i)-y2_dom(j))], inf)*norm(shape, inf)* ...
                        norm([-2*(x1_dom(i)-y1_dom(j)) 0; 0 -2*(x2_dom(i)-y2_dom(j))], inf) ...
                        + 2*norm((X_dom(2,i)-Y_dom(2,j).^2)', inf)*norm(shape, inf)*norm([-2 0; 0 -2], inf), ...
                        2*norm((X_dom(2,i)-Y_dom(2,j).^2)', inf)*norm(shape, inf)*norm([-2*(x1_dom(i)-y1_dom(j)) 0; 0 -2*(x2_dom(i)-y2_dom(j))], inf) + ...
                        2*norm((X_dom(2,i)-Y_dom(2,j).^2)', inf)*norm(shape, inf)*norm([-2*(x1_dom(i)-y1_dom(j)) 0; 0 -2*(x2_dom(i)-y2_dom(j))], inf)];
            dE2_dxdy = [dE2_dxdy, temp];
        end
    end

    A = norm(dh_dx, inf); %*0.01;
    B = norm(dE2_dxdy, inf)*norm([1 0; 1 0; 0 1], inf); %*0.1;
    E = [];
    U = [];


    controller = @(x, vel_des) Sampling_Controller(x, vel_des, B, A, shape, P_safe, gamma, num_samp, tau);
end
%% Setup for the SOS controller
if strcmp(cont, 'sos')
    syms y1 y2 u1 u2 t real;
    syms theta x1 x2 real;
    y = [y1; y2]; u_sim = [u1; u2];
    del = ([x1; x2] - y);
    rot = [cos(theta), sin(theta) ; -sin(theta), cos(theta)];
    new_coords = shape*rot*del;                                           
    Vsym = sum(new_coords.^4) - (shape(1,1)*shape(2,2))^4;                           
    diffVsym = [diff(Vsym,x1), diff(Vsym,x2), diff(Vsym, theta)];
    controller = @(x, vel_des) SOS_controller(x, vel_des, Vsym, diffVsym, x1, x2, theta, y, u_sim, t, P_safe); 
end

%% Setup for the point controller
if strcmp(cont, 'point')
    gamma = 0.6;
    controller = @(x, vel_des) Point_controller(x, vel_des, P_safe_r, gamma);
end

%% Start the robotarium simulation
% Initialize the robotarium
r = Robotarium('NumberOfRobots', N, 'ShowFigure', true);
automatic_parking_controller = create_automatic_parking_controller();

% Plot the safe Set
plot_safeSet(P_safe, 'b')
hold on
plot_safeSet(P_safe_r, 'r')

for t = 0:1000
   x = r.get_poses();
   dxi(:, 1) = [0 ; 0];
   dxi = automatic_parking_controller(x, [0.2; -0.3; 0]);
   r.set_velocities(1:N, dxi);
   r.step();
end

x = r.get_poses();
Plt_data1 = [];
Plt_data1 = [Plt_data1; x(1,1); x(2,1)];
p1 = plot(Plt_data1(1), Plt_data1(2), 'k-.', 'LineWidth', 3);
drawnow

vol = plot_squircle(x, shape);
drawnow
grid on
r.step();
T_comp = 0;
T = [];
t_sim = 0;
U = [];

for t = 0:iterations
    delete(vol);
    x = r.get_poses();
    vel_des = [1; 0.4];
    
    if strcmp(cont, 'extent')
        % Sampling based controller 
        [dxu, comp_time, ext] = controller(x, vel_des);
        E = [E, ext];
        dxu = 0.1*dxu;
    elseif strcmp(cont, 'sos')
        % SOS controller 
        [dxu, comp_time] = controller(x, vel_des);
        dxu = 1e5*double(dxu);
    elseif strcmp(cont, 'point')
        % Single point controller
        [dxu, comp_time] = controller(x, vel_des);
        dxu = 0.25*dxu;
    elseif strcmp(cont, 'none')
        dxu = vel_des;
        comp_time= 0;
    end
    
    
    T_comp = T_comp + comp_time;
    U = [U, norm(dxu)];
    
    % Scale the output
    r.set_velocities(1:N, dxu);
    r.step();
    
    % Update plot of trajectory and Extent set
    Plt_data1 = [Plt_data1, [x(1,1); x(2,1)]];
    p1.XData = Plt_data1(1,:);
    p1.YData = Plt_data1(2,:);
    t_sim = t_sim + (1/30);
    T = [T, t_sim];
    vol = plot_squircle(x, shape);
    drawnow 
end


r.debug();