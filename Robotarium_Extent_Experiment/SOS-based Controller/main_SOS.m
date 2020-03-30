% Extent-Compatible Control Barrier Functions
% NOTE: Please run the "init.m" file from the Robotarium package prior to
%       executing this code
% Authors: Mohit Srinivasan, Matt Abate, Gustav Nilsson, and Samuel Coogan
% Date Modified: 03/25/2020

clear all;
close all;
clc;

N = 1;
iterations = 1000;
r = Robotarium('NumberOfRobots', N, 'ShowFigure', true);
automatic_parking_controller = create_automatic_parking_controller();

U = [];
global P_safe
P_safe = [1/1^2 0; 0 1/0.8^2];

global shape
shape = [1 0; 0 1.3333];

a = 0; b = 0;
plot_safeSet(P_safe, a, b, 'b');
hold on

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

r.step();

%% SOS Setup
syms y1 y2 u1 u2 t real;
syms theta x1 x2 real;
y = [y1; y2]; u_sim = [u1; u2];
del = ([x1; x2] - y);
rot = [cos(theta), sin(theta) ; -sin(theta), cos(theta)];
new_coords = shape*rot*del;                                           
Vsym = sum(new_coords.^4) - (0.1333)^4;                               
diffVsym = [diff(Vsym,x1), diff(Vsym,x2), diff(Vsym, theta)];
    
vol = plot_squircle(x);
drawnow
hold on
grid on
T_comp = 0;
T = [];
t_sim = 0;

for i = 0:iterations
    
    delete(vol);
    x = r.get_poses();
    vel_des = [1; 0.4];
    
    % Extent-Compatible Barrier Functions Algorithm
    [dxu, comp_time] = SOS_controller(x, vel_des, Vsym, diffVsym, x1, x2, theta, y, u_sim, t);
    t_sim = t_sim + comp_time;
    dxu = 1e5*double(dxu);
    U = [U, norm(dxu)];
    r.set_velocities(1:N, dxu);
    r.step();
    
    % Update plot of trajectory and Extent set
    Plt_data1 = [Plt_data1, [x(1,1); x(2,1)]];
    p1.XData = Plt_data1(1,:);
    p1.YData = Plt_data1(2,:);

    vol = plot_squircle(x);
    drawnow
        
end

display(t_sim/iterations);
r.debug();