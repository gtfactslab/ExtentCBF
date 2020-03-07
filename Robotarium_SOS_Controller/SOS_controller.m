function y = SOS_controller(x, u_nom)

    % Creat symbolic variables for SOS implementation
    syms y1 y2 u1 u2 t real;
    syms theta x1 x2 real;
    y = [y1; y2]; u_sim = [u1; u2];
    del = ([x1; x2] - y);
    rot = [cos(theta), sin(theta) ; -sin(theta), cos(theta)]; % Rotation matrix
    new_coords = [1.5 0; 0 2]*rot*del;
    Vsym = sum(new_coords.^4) - (0.2)^4; % Extent function
    diffVsym = [diff(Vsym,x1), diff(Vsym,x2), diff(Vsym, theta)];
    P_safe = [1/0.8^2 0; 0 1/0.8^2];
    
    % Setup SOS
    sostoolsoptions.solver = 'sdpt3';
    sostoolsoptions.frlib.approx = 'dd';
    sostoolsoptions.frlib.useQR = 1;

    % Parameters for safe set, Extent Set and sampling controller
    vd = u_nom;    
    V = subs(Vsym,{x1, x2, theta},{x(1),x(2), x(3)});
    dVdx = subs(diffVsym,{x1, x2, theta},{x(1),x(2), x(3)});
    LgV = dVdx * [cos(x(3)), 0; sin(x(3)), 0; 0, 1];
    LgV = vpa(LgV);
                
    monos = monomials(y, [0 1 2]); %Generate vector of monomials
    h_of_y = 1 - y'*P_safe*y;
    
    %Run SOS proram
    u = u_sim;
    prog = sosprogram(y);
    prog = sosdecvar(prog,u);
    prog = sosdecvar(prog,t);
    [prog,s2] = sossosvar(prog,monos);
    [prog,s3] = sossosvar(prog,monos);
    prog = sosineq(prog, LgV*u + 150*h_of_y + 150.*V - s2*h_of_y - s3*V);
    prog = sossetobj(prog,t);
    M = [eye(2), u; u',t-vd'*vd+2*vd'*u];
    prog = sosmatrixineq(prog,M,'quadraticMineq');
    prog = sossolve(prog,sostoolsoptions);
    y = sosgetsol(prog,u);
        
end