function [LgE, min_E, y, comp_time] = Sampling_Controller(x, u_nom, L_Edot, L_h)
    
    % Parameters for the safe set, extent Set and sampling controller
    gamma = 0.8;
    tau = pi/1000;
    P_safe = [1/1^2 0; 0 1/0.8^2];
    c1 = 0; c2 = 0;
    C = [c1; c2];
    
    shape = [1.5, 0; 0, 2];
    th = 0:tau:2*pi;
    H = eye(2);                              
    a = [];
    b = [];
    E_sys = [];
    
    phi = x(3);
    rot = [cos(phi) sin(phi); -sin(phi) cos(phi)]*[0 1; -1 0];

    % Sampling based constraints       
    for i=1:length(th)
        
        test_units = [sqrt(abs(cos(th(i))))*0.15*sign(cos(th(i))); ...
                        sqrt(abs(sin(th(i))))*0.2*sign(sin(th(i)))];
                    
        testy = (test_units'*rot)' + x(1:2);
        del = (x(1:2) - testy).^2;
        
        dEdx = [2*del'*rot'*shape*rot*2*[(x(1)-testy(1)) 0; 0 (x(2)-testy(2))], ...
                del'*[-sin(phi), -cos(phi); cos(phi), -sin(phi)]*shape*rot*del + ...
                del'*rot*shape*[-sin(phi), -cos(phi); cos(phi), -sin(phi)]*del];
        LgE = dEdx * [cos(phi), 0; sin(phi), 0; 0, 1];
        
        h_safe_of_y = 1 - (testy - C)'*P_safe*(testy - C);
        E_sys = [E_sys, h_safe_of_y];
        a = [a; -LgE];
        b = [b; gamma*(h_safe_of_y) - (gamma*L_h + L_Edot)*tau];

    end
    
    % Solve QP
    min_E = min(E_sys);
    opts = optimoptions(@quadprog, 'Display', 'off', 'TolFun', 1e-5, 'TolCon', 1e-4);
    tic
    y = quadprog(sparse(H), double(-2*u_nom'), a, b, [], [], [], [], [], opts);  
    comp_time = toc;
end
