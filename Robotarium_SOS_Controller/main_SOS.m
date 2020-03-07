% Extent Compatible Control Barrier Functions for Continuous Time Dynamical
% Systems
% NOTE: Please run the "init.m" file from the Robotarium package prior to
%       executing this code
% Authors: Mohit Srinivasan, Matt Abate, Gustav Nilsson, and Samuel Coogan
% 03/06/2020

clear all;
close all;
clc;

N = 1;
iterations = 2000;
r = Robotarium('NumberOfRobots', N, 'ShowFigure', true);

automatic_parking_controller = create_automatic_parking_controller();
uni_barrier_certificate = create_uni_barrier_certificate2();

U = [];
                 
P = [1/0.8^2 0; 0 1/0.8^2];
a = 0; b = 0;
plot_safeSet(P, a, b, 'r');
hold on

for t = 0:1000
    
   x = r.get_poses();
   dxi(:, 1) = [0 ; 0];
   dxi = automatic_parking_controller(x, [0.5; 0; 0]);
   [dxi, ~] = uni_barrier_certificate(dxi, x, []);
   r.set_velocities(1:N, dxi);
   r.step();
    
end

x = r.get_poses();
Plt_data1 = [];
Plt_data1 = [Plt_data1; x(1,1); x(2,1)];
p1 = plot(Plt_data1(1), Plt_data1(2), 'k-.', 'LineWidth', 3);
drawnow

syms y1 y2 u1 u2 t real;
syms theta x1 x2 real;
y = [y1; y2]; u_sim = [u1; u2];
del = ([x1; x2] - y);
rot = [cos(theta), sin(theta) ; -sin(theta), cos(theta)]; % Rotation matrix
new_coords = [1.5 0; 0 2]*rot*del;
Vsym = sum(new_coords.^4) - (0.2)^4; % Extent function

P = [1/1.5^2 0; 0 1/1^2]/0.1^2;
vol = plot_squircle(x, Vsym, x1, x2, theta);
drawnow
grid on
r.step();

for t = 0:iterations
    
    delete(vol);
    x = r.get_poses();
    vel_des = [1; 1];
    
    % Extent-Compatible Barrier Functions Algorithm
    dxu = SOS_controller(x, vel_des);
    U = [U, norm(dxu)];
    dxu = uni_barrier_certificate(dxu, x);
    dxu = 1*dxu;
    r.set_velocities(1:N, dxu);
    r.step();
    
    % Update plot of trajectory and Extent set
    Plt_data1 = [Plt_data1, [x(1,1); x(2,1)]];
    p1.XData = Plt_data1(1,:);
    p1.YData = Plt_data1(2,:);

    vol = plot_squircle(x, Vsym, x1, x2, theta);
    drawnow
    
        
end

r.debug();