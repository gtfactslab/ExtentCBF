function [LgE, min_E, y] = extent_barr_ellipse_Lp(x, u_nom, L_Edot, L_h)
    
    % Parameters for the Lp safe set, Extent Set and sampling controller
    gamma = 5;
    tau = pi/100;
    p = 4;
    sigma = [1, 0.4, pi/20];
    theta_k = pi/2;
    k = theta_k/(2*sigma(1));
    c = 1;
    theta0 = sign(k)*pi/2;
    
    shape = [1.5, 0; 0, 1];
    
    th = 0:tau:2*pi;
    test_units = [0.15*cos(th); 0.1*sin(th)];
    
    H = eye(2);
    
    % Weighted Polar Lp Safe Set
    R = sqrt((k*x(1))^2 + (k*x(2) + 1)^2);
    theta = atan2((k*x(2) + 1), k*x(1));
    
    alpha = (R - c)/sigma(2);
    beta = (theta - theta0)/sigma(1);
    
    paralpha = [(k^2*x(1))/(sigma(2)*R), ...
               (k*(k*x(2)+1))/(sigma(2)*R)];
          
    parbeta = [-(k*(k*x(2)+1))/(sigma(1)*R^2), ....
               (k^2*x(1))/(sigma(1)*R^2)];   
                                    
    a = [];
    b = [];
    E_sys = [];
    
    phi = x(3);
    rot = [cos(phi) sin(phi); -sin(phi) cos(phi)];
    drot = [-sin(phi) cos(phi); -cos(phi) -sin(phi)];
    
    P=rot'*shape*rot;
    T = x(3);
    
    % Sampling based constraints       
    for i=1:length(th)
        
        testy = x(1:2) + P*test_units(:,i);
        del = x(1:2) - testy;
        
        dE_dx = [2*del'*rot'*shape*rot, del'*(drot'*shape*rot+rot'*shape*drot)*del];
         
        LgE = dE_dx*[cos(phi) 0; sin(phi) 0; 0 1];
                
        R_of_y = sqrt((k*testy(1))^2 + (k*testy(2) + 1)^2);
        theta_of_y = atan2((k*testy(2) + 1), k*testy(1));

        alpha_of_y = (R_of_y - c)/sigma(2);
        beta_of_y = (theta_of_y - theta0)/sigma(1);

        h_safe_of_y = abs(k) - (alpha_of_y^p + beta_of_y^p)^(1/p);                                          
        
        E_sys = [E_sys, h_safe_of_y];
        a = [a; -LgE];
        b = [b; gamma*(h_safe_of_y) - (gamma*L_h + L_Edot)*tau];

    end
    
    % Solve QP
    min_E = min(E_sys);
    opts = optimoptions(@quadprog, 'Display', 'off', 'TolFun', 1e-5, 'TolCon', 1e-4);
    y = quadprog(sparse(H), double(-2*u_nom'), a, b, [], [], [], [], [], opts);  
        
end