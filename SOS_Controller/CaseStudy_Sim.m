% Extent Compatible Control Barrier Function Paper - Case Study
% Unicycle Surrounded by Squircle
% Author: Mohit Srinivasan, Matthew Abate, Gustav Nilsson, Samuel Coogan
% Date: 12/23/2019

% this script was written for MATLAB 2019b, and requires SOSTOOLS with
% solver SDPT3

%addsostools
%installmex(1)
%startup

close all; 
clear all;
clc; 

% ------------------------------------------------
% Generate control inputs and initial conditions
% ------------------------------------------------

tstep=.02; time=0:tstep:2;
n=length(time);             % number of simulated time steps
options = sdpsettings('verbose',0,'solver','sdpt3');

vel_des = ones(n, 1);       % Constant desired velocity
rot_des = ones(n, 1);       % Constant desired rotation

u_des = [vel_des, rot_des]; % Concatinate desired input
clear vel_des rot_des


% The initial states are rotationally symeteric about (x, y) = (0, 0)
dist = .5;        % distance of unicycles to origin
angle_offset = 0; % initial angle offset for each unicycle

% Get initial conditions where:
%  - NBF means no barrier function
%  - BF means traditional barrier function
%  - SOS means extent compatible barrier function
x_NBF = [dist, 0, pi/2 + angle_offset];
x_BF  = [dist*cos(2*pi/3), dist*sin(2*pi/3), 7*pi/6 + angle_offset];
x_SOS = [dist*cos(4*pi/3), dist*sin(4*pi/3), 11*pi/6 + angle_offset];


% Set-up figure
figure(1); clf;
hold on; grid on;

circle(0,0,1); % plot safe set
axis([-1.5, 1.5, -1.5, 1.5])
pbaspect([1 1 1])
xlabel('$x_1$','Interpreter','latex')
xticks([-1, 0, 1])
ylabel('$x_2$','Interpreter','latex')
yticks([-1, 0, 1])
set(gca,'TickLabelInterpreter','latex', 'FontSize', 16)
legend('Location', 'southoutside')
drawnow



% ------------------------------------------------
% Setup Simulation
% ------------------------------------------------

% Setup SOS
sostoolsoptions.solver = 'sdpt3';
sostoolsoptions.frlib.approx = 'dd';
sostoolsoptions.frlib.useQR = 1;

% Creat symbolic variables for SOS implementation
syms y1 y2 u1 u2 t real;
syms theta x1 x2 real;
y = [y1; y2]; u_sim=[u1;u2];
xs = [x1; x2; theta];
del = ([x1; x2] - y);
rot = [cos(theta), sin(theta) ; -sin(theta), cos(theta)]; % Rotation matrix
new_coords = [1.5 0; 0 2]*rot*del;
Vsym = sum(new_coords.^4) - (0.2)^4; % Extent function

% Create object for storing plot objects
plot_holder = [];

% ------------------------------------------------
% RUN SIMULATION
% ------------------------------------------------
for i = 1:n
    hold on
    vd = u_des(i,:)';      % current suggested input
        
    % ------------------------------------------------
    % Simulate unicycle dynamics with no barrier
    % ------------------------------------------------
    x_curr = x_NBF(end,:);
    x_next = [];
    x_next(1, 1) = x_curr(1, 1) + tstep * vd(1) * cos(x_curr(1, 3));
    x_next(1, 2) = x_curr(1, 2) + tstep * vd(1) * sin(x_curr(1, 3));
    x_next(1, 3) = x_curr(1, 3) + tstep * vd(2);
    x_NBF=[x_NBF;x_next(1, :)]; %Store
    
    
    % ------------------------------------------------
    % Simulate dynamics with traditional barrier
    % ------------------------------------------------
    % NOTE: In this case, we model the unicycle as single integrator, using the
    %       tradidional diffeomorphism. 
    x_curr = x_BF(end,:)'; % current state
    T = x_curr(3);          % current angle state
    
    % Diffeomorphism
    l = .05; % parameter
    A = [cos(T), -l*sin(T); sin(T), l*cos(T)];
    z =[x_curr(1)+l*cos(x_curr(3)) ; x_curr(2)+l*sin(x_curr(3))];
    
    % Compute Barrier
    Lgh = -2*z' * A; % Safe set: h = (1 - z'*z)
    u = sdpvar(2,1);
    Constraints = [-Lgh*u <= 5*(1 - z'*z)];
    optimize(Constraints,(u-vd)'*(u-vd), options);
    u_act = value(u);
    
    % Simulate
    x_next = [];
    x_next(1, 1) = x_curr(1) + tstep * u_act(1) * cos(x_curr(3, 1));
    x_next(2, 1) = x_curr(2) + tstep * u_act(1) * sin(x_curr(3, 1));
    x_next(3, 1) = x_curr(3) + tstep * u_act(2);
    x_BF = [x_BF; x_next'];
    
    
    % ------------------------------------------------
    % Simulate dynamics with extent-compatible barrier
    % ------------------------------------------------
    x_curr = x_SOS(end,:)';                  % Current State
    T = x_curr(3);                            % Current angle state
    h_eval = 1 - x_curr(1:2)' * x_curr(1:2);  % Current Barrier evaluation 
    
    %Symbolic representation of super ellipsoidal extent function    
    V = subs(Vsym,{x1, x2, theta},{x_curr(1),x_curr(2), x_curr(3)});
    diffVsym = [diff(Vsym,x1), diff(Vsym,x2), diff(Vsym, theta)];
    dVdx = subs(diffVsym,{x1, x2, theta},{x_curr(1),x_curr(2), x_curr(3)});
    LgV = dVdx * [cos(T), 0; sin(T), 0; 0, 1];
    LgV = vpa(LgV);
                
    monos = monomials([y], [0 1 2]); %Generate vector of monomials
    h_of_y = 1-y'*y;
    
    %Run SOS proram
    u = u_sim;
    prog = sosprogram([y1,y2]);
    prog = sosdecvar(prog,u);
    prog = sosdecvar(prog,t);
    [prog,s2] = sossosvar(prog,monos);
    [prog,s3] = sossosvar(prog,monos);
    prog = sosineq(prog, LgV*u + 150*h_of_y + 150.*V - s2*h_of_y-s3*V);
    prog = sossetobj(prog,t);
    M = [eye(2), u; u',t-vd'*vd+2*vd'*u];
    prog = sosmatrixineq(prog,M,'quadraticMineq');
    prog = sossolve(prog,sostoolsoptions);
    u_act = sosgetsol(prog,u);

    % Simulate Unicycle Dynamics
    x_next = [];
    x_next(1, 1) = x_curr(1) + tstep * u_act(1) * cos(x_curr(3, 1));
    x_next(2, 1) = x_curr(2) + tstep * u_act(1) * sin(x_curr(3, 1));
    x_next(3, 1) = x_curr(3) + tstep * u_act(2);
    x_SOS=[x_SOS; x_next'];
    
    % Check to make sure SOSTOOLS did not return an error
    if prog.solinfo.info.pinf == 1 || prog.solinfo.info.dinf == 1
        error('SOStools did not find a solution\n')
    end

    % Update plot
    
    colors = [1 0 0; 0 1 0; 0 0 1];
    
    % plot unicycle with no barrier
    delete(plot_holder)
    plot_holder(1, 1) = plot(x_NBF(:,1),x_NBF(:,2), ...
                            'Color', .5*colors(1, :), ...
                            'LineWidth', 1.5, ...
                            'DisplayName', 'Nominal trajectory');
                        
    % plot unicycle with traditional barrier
    plot_holder(1, 2) = plot(x_BF(:,1),x_BF(:,2), ...
                            'Color', .5*colors(3, :), ...
                            'LineWidth', 1.5, ...
                            'DisplayName', 'ZCBF trajectory');
         
    % plot unicycle with extent barrier  
    plot_holder(1, 3) = plot(x_SOS(:,1),x_SOS(:,2), ...
                            'Color', .5*colors(2, :), ...
                            'LineWidth', 1.5, ...
                            'DisplayName', 'Ex-CBF trajectory');
      
    plot_holder(2:3, 1) = plot_ellipse(x_NBF(end,:), colors(1, :));
    plot_holder(2:3, 2) = plot_ellipse(x_BF(end,:), colors(3, :));
    plot_holder(2:3, 3) = plot_ellipse(x_SOS(end, :), colors(2, :));
    
    drawnow  % plot
    time(i)  % display the time
end
clear x_next x_curr vd l t A z



% This function plots an ellipse
% input: unicycle state (3x1 vector) and color (1x3 vector)
% output: plot handle
function out = plot_ellipse(x_curr, col)
    syms y1 y2 u1 u2 t real;
    syms theta x1 x2 real;
    y=[y1; y2]; u=[u1;u2];
    xs=[x1; x2; theta];
    
    %Super Ellipsoidal extent function
    rot = [cos(theta), sin(theta) ; -sin(theta), cos(theta)]; % Rotation matrix
    del = ([x1; x2] - y);
    new_coords=[1.5 0; 0 2]*rot*del;
    Vsym = sum(new_coords.^4) - (0.2)^4;
    
    V=subs(Vsym,{x1, x2, theta}, {x_curr(1), x_curr(2), x_curr(3)});
    
    fc=fcontour(matlabFunction(V),[-1.2,1.2], 'HandleVisibility', 'off');
    fc.LevelList=[0];
    fc.LineColor = 'r';
    fc.LineWidth = 1.2;
    drawnow

    [cc, hh] = contour(get(fc,'XData'), get(fc,'YData'), get(fc,'ZData'), ...
                            [0 0], 'r', ...
                            'HandleVisibility', 'off'); % Overlay contour line
    delete(fc)
    fc = patch(cc(1,2:end), cc(2, 2:end), col, ...
                    'EdgeColor', .5*col, ...
                    'FaceAlpha', .3, ...
                    'HandleVisibility', 'off');
    pause(.01);
    delete(hh);
    
    % Plot pint triangle
    r = .1;     
    q = [x_curr(1, 1) + r/1.5 * cos(x_curr(1, 3)), ...
         x_curr(1, 2) + r/1.5 * sin(x_curr(1, 3)); ...
         x_curr(1, 1) + r/3 * cos(x_curr(1, 3) + 2*pi/3), ...
         x_curr(1, 2) + r/3 * sin(x_curr(1, 3) + 2*pi/3); ...
         x_curr(1, 1) + r/3 * cos(x_curr(1, 3) + 4*pi/3), ...
         x_curr(1, 2) + r/3 * sin(x_curr(1, 3) + 4*pi/3)];
    qq = patch(q(:, 1), q(:, 2), 'm', 'HandleVisibility', 'off');
    
    out = [fc; qq];
end
