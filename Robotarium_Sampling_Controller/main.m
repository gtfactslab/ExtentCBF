% Extent Compatible Control Barrier Functions for Continuous Time Dynamical
% Systems
% NOTE: Please run the "init.m" file from the Robotarium package prior to
%       executing this code
% Authors: Mohit Srinivasan, Matt Abate, Gustav Nilsson, and Samuel Coogan
% 01/20/2020

clear all;
close all;
clc;

N = 1;
iterations = 2000;
r = Robotarium('NumberOfRobots', N, 'ShowFigure', true);

automatic_parking_controller = create_automatic_parking_controller();
uni_barrier_certificate = create_uni_barrier_certificate2();

L_Edot = norm([1 0; 1 0; 0 1], inf)*0.1;
L_h = 4;
E = [];
U = [];
                          
% Plot the Lp Safe Set
P_safe = [1/0.8^2 0; 0 1/0.8^2];
c1 = 0; c2 = 0;
plot_safeSet(P_safe, c1, c2, 'r')
hold on

for t = 0:1000
    
   x = r.get_poses();
   dxi(:, 1) = [0 ; 0];
   dxi = automatic_parking_controller(x, [-0.2; 0; 0]);
   [dxi, ~] = uni_barrier_certificate(dxi, x, []);
   r.set_velocities(1:N, dxi);
   r.step();
    
end

x = r.get_poses();
Plt_data1 = [];
Plt_data1 = [Plt_data1; x(1,1); x(2,1)];
p1 = plot(Plt_data1(1), Plt_data1(2), 'k-.', 'LineWidth', 3);
drawnow

P = [1/1.5^2 0; 0 1/1^2]/0.1^2;
vol = plot_squircle(x);
drawnow
grid on
r.step();

for t = 0:iterations
    
    delete(vol);
    x = r.get_poses();
    vel_des = [0.1; 0.2];
    
    % Extent-Compatible Barrier Functions Algorithm
    [LgE, ext, dxu] = Sampling_Controller(x, vel_des, L_Edot, L_h);
    E = [E, ext];
    U = [U, norm(dxu)];
    dxu = uni_barrier_certificate(dxu, x);
    r.set_velocities(1:N, dxu);
    r.step();
    
    % Update plot of trajectory and Extent set
    Plt_data1 = [Plt_data1, [x(1,1); x(2,1)]];
    p1.XData = Plt_data1(1,:);
    p1.YData = Plt_data1(2,:);

    vol = plot_squircle(x);
    drawnow
    
        
end

r.debug();