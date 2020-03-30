function [y, comp_time] = SOS_controller(x, u_nom, Vsym, diffVsym, x1, x2, theta, y, u_sim, t)
    
    global P_safe
    
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
    prog = sosineq(prog, LgV*u + 0.1*h_of_y^6 + 0.1.*V - s2*h_of_y - s3*V);
    prog = sossetobj(prog,t);
    M = [eye(2), u; u',t-vd'*vd+2*vd'*u];
    prog = sosmatrixineq(prog,M,'quadraticMineq');
    tic
    prog = sossolve(prog,sostoolsoptions);
    y = sosgetsol(prog,u);
    comp_time = toc;
    
end