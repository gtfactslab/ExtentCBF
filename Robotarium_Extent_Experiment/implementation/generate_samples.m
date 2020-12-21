function theta = generate_samples(shape,tau)  
    

    options =  optimset('TolFun',1e-12, 'Display', 'off');

    initial_stepsize = tau/max(diag(shape));
    
    theta = [0];
    last_theta = 0;
    stop = 0;
    while stop == 0
        samp_dist = @(x) (sqrt(abs(cos(last_theta))).*shape(1, 1).*sign(cos(last_theta))- sqrt(abs(cos(x))).*shape(1, 1).*sign(cos(x))).^2 + (sqrt(abs(sin(last_theta))).*shape(2, 2).*sign(sin(last_theta)) - sqrt(abs(sin(x))).*shape(2, 2).*sign(sin(x))).^2 - tau^2;
        i = 2;
        new_theta = fsolve(samp_dist, last_theta+i*initial_stepsize,options);
        
        while new_theta < last_theta + 1e-8
            i = i*2;
            new_theta = fsolve(samp_dist, last_theta+i*initial_stepsize,options);
        end

        if new_theta > pi/2 %new_theta > last_theta+4*initial_stepsize ||
            stop = 1;
        else
            last_theta = new_theta;
            theta = [theta; last_theta];
        end
    end
    theta = [theta; flip(pi-theta); theta+pi; flip(2*pi-theta)];
end


