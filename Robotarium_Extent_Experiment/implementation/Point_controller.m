function [y, comp_time, min_E] = Point_controller(x, u_nom, P_safe, gamma, shape, N)
    xpos = [x(1); x(2)];
    phi = x(3);

    H = eye(2);
    g = [cos(phi) 0; sin(phi) 0; 0 1];
    Lg = [-2*xpos'*P_safe 0];
    h = 1-xpos'*P_safe*xpos;

    % Set up QP
    a = -Lg*g;
    b = gamma * h;

    % Solve QP
    opts = optimoptions(@quadprog, 'Display', 'off', 'TolFun', 1e-5, 'TolCon', 1e-4);
    tic
    y = quadprog(sparse(H), double(-2*u_nom'), a, b, [], [], [], [], [], opts);
    comp_time = toc;
    
    
  
    th = 0:(2*pi/N):2*pi;
    phi = x(3);
    rot = [-sin(phi) cos(phi); -cos(phi) -sin(phi)];
    E_sys = [];
    for i=1:length(th)
        
        test_units = [sqrt(abs(cos(th(i))))*shape(1, 1)*sign(cos(th(i))); ...
            sqrt(abs(sin(th(i))))*shape(2, 2)*sign(sin(th(i)))];
        
        testy = (test_units'*rot)' + x(1:2);
       
        h_safe_of_y = 1 - (testy)'*P_safe*(testy);
        E_sys = [E_sys, h_safe_of_y];
    end
    
    % Solve QP
    min_E = min(E_sys);
end