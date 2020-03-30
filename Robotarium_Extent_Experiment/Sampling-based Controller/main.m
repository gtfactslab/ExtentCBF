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

% Plot the safe Set
global P_safe
P_safe = [1/1^2 0; 0 1/0.8^2];
c1 = 0; c2 = 0;
plot_safeSet(P_safe, c1, c2, 'b')
hold on

x1_dom = -1.6:0.01:1.6; x2_dom = -1.6:0.01:1.6;
y1_dom = -1.6:0.01:1.6; y2_dom = -1.6:0.01:1.6;

global shape
shape = [1 0; 0 1.333];

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

A = norm(dh_dx, inf)*0.001;
B = norm(dE2_dxdy, inf)*norm([1 0; 1 0; 0 1], inf)*0.01;
E = [];
U = [];

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

vol = plot_squircle(x);
drawnow
grid on
r.step();
T_comp = 0;
T = [];
t_sim = 0;

for t = 0:iterations
    
    delete(vol);
    x = r.get_poses();
    vel_des = [1; 0.4];
    
    % Extent-Compatible Barrier Functions Algorithm
    [LgE, ext, dxu, comp_time] = Sampling_Controller(x, vel_des, B, A);
    T_comp = T_comp + comp_time;
    E = [E, ext];
    U = [U, norm(dxu)];
    dxu = 0.25*dxu;
    r.set_velocities(1:N, dxu);
    r.step();
    
    % Update plot of trajectory and Extent set
    Plt_data1 = [Plt_data1, [x(1,1); x(2,1)]];
    p1.XData = Plt_data1(1,:);
    p1.YData = Plt_data1(2,:);
    t_sim = t_sim + (1/30);
    T = [T, t_sim];
    vol = plot_squircle(x);
    drawnow
  
end

% display(T_comp/iterations);
% display(min(E));
% save('data_50.mat', 'E', 'T');
% save('comp_time.mat', 'T_comp')

r.debug();