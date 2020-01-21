function d = checkDonut(x, y, nat)

p = 4;
sigma = [1, 0.4];
theta_k = pi/2;
k = theta_k/(2*sigma(1));
c = 1;
theta0 = sign(k)*pi/2;

R = ((k*x)^2 + (k*y+1)^2)^(1/2);
theta = atan2(k*y+1, k*x);

alpha = abs(R - c)/sigma(2);
beta = abs(theta - theta0)/sigma(1);

Hg = abs(k) - (alpha^p + beta^p)^(1/p);

if (nat == 1)
    if(Hg >= 0)
        d = 1;
    else
        d = 0;
    end
else
    if(Hg >= 0)
        d = 0;
    else
        d = 1;
    end

end