function [LgE, min_E, y, comp_time] = Sampling_Controller(x, u_nom, B, A, shape, P_safe, gamma, N)
        
    % Parameters for the safe set, extent set and sampling controller 
   
    c1 = 0; c2 = 0;
    C = [c1; c2];
    
    tau = (2*pi/N);
    th = 0:tau:2*pi;
    H = eye(2);
    a = []; b = [];
    E_sys = [];
    
    phi = x(3);
    rot = [-sin(phi) cos(phi); -cos(phi) -sin(phi)];
    drot = [-cos(phi) -sin(phi); sin(phi) -cos(phi)];
    
    rot_trans = [-sin(phi) -cos(phi); cos(phi) -sin(phi)];
    drot_trans = [-cos(phi) sin(phi); -sin(phi) -cos(phi)];
    
    % Sampling based constraints
    for i=1:length(th)
        
        test_units = [sqrt(abs(cos(th(i))))*shape(1, 1)*sign(cos(th(i))); ...
                        sqrt(abs(sin(th(i))))*shape(2, 2)*sign(sin(th(i)))];
                                        
                    
        testy = (test_units'*rot)' + x(1:2);
        del = (x(1:2) - testy).^2;
        
        dEdx = [2*del'*rot_trans*shape*rot*2*[(x(1)-testy(1)) 0; 0 (x(2)-testy(2))], ...
                del'*drot_trans*shape*rot*del + ...
                del'*rot'*shape*drot*del];
        
        LgE = dEdx * [cos(phi), 0; sin(phi), 0; 0, 1];
        
        h_safe_of_y = 1 - (testy - C)'*P_safe*(testy - C);
        E_sys = [E_sys, h_safe_of_y];
        a = [a; -LgE];
        b = [b; gamma*(h_safe_of_y) - (gamma*A + B)*tau];
        
    end
    
    % Solve QP
    min_E = min(E_sys);
    opts = optimoptions(@quadprog, 'Display', 'off', 'TolFun', 1e-5, 'TolCon', 1e-4);
    tic
    y = quadprog(sparse(H), double(-2*u_nom'), a, b, [], [], [], [], [], opts);  
    comp_time = toc;
    
end