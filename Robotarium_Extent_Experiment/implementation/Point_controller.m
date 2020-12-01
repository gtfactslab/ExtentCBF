function [y, comp_time] = Point_controller(x, u_nom, P_safe, gamma)
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
end