function max_grad = findMinGrad_Donut()

p = 4;
sigma = [1, 0.4, pi/20];
theta_k = pi/2;
k = theta_k/(2*sigma(1));
c = 1;
theta0 = sign(k)*pi/2;
epsilon = 1;

x_Lp = linspace(-1.6, 1.6, 100);
y_Lp = linspace(-1, 1, 100);

[x, y] = meshgrid(x_Lp, y_Lp);

x_new = k.*x;
y_new = k.*y + 1;

R = sqrt(((x_new).^2 + (y_new).^2));
R_k = 1/abs(k);
theta = atan2(y_new, x_new);
alpha_b = sign(k);

Bx = (R_k + y).*cos(alpha_b*pi/2 + (theta_k/2).*(x./sigma(1)));
By = (R_k + y).*sin(alpha_b.*pi/2 + (theta_k/2).*(x./sigma(1))) - alpha_b.*R_k;

D = [];
Dx = [];
Dy = [];
m = 1;
for i = 1:length(Bx(:,1))
    
    for j = 1:length(Bx(1,:))
        donut = checkDonut(Bx(i,j), By(i,j), 1);
        if(donut == 1)
            Dx = [Dx, Bx(i,j)];
            Dy = [Dy, By(i,j)];
            m = m+1;
            %plot(Bx(i,j), By(i,j), 'bo');
            %hold on
        end
    end
end

D = [Dx', Dy'];
G = zeros(length(Dx), 2);
Norm_Grad = [];

for i = 1:length(D(:,1))
    
    x_new_B = k*D(i,1);
    y_new_B = k*D(i,2)+1;
    R_B = sqrt(((x_new_B)^2 + (y_new_B)^2));
    theta_B = atan2(y_new_B, x_new_B);
    Lp_weight_B_p_1 = 0.5*((R_B - c)^(p)/(sigma(2)^p)+(theta_B - theta0)^(p)/(sigma(1)^p))^(1/p-1);
    
    V_v = Lp_weight_B_p_1*p*[(R_B - c)^(p-1)/(sigma(2)^p); -(theta_B - theta0)^(p-1)/(sigma(1)^p)];
    Q_inv = [cos(theta_B), sin(theta_B); (1/R_B)*sin(theta_B), -(1/R_B).*cos(theta_B)];
    G(i, :) = V_v'*Q_inv; 
    Norm_Grad = [Norm_Grad, norm(G(i,:))];
    
end

% plot3(D(:,1)', D(:,2)', Norm_Grad, 'ro');
xlim([-1.6 1.6])
ylim([-1 1])

max_grad = max(Norm_Grad);

end
