function plot_safeSet(P, a, b, c)

      theta = 0:0.001:2*pi;
      x = (1/sqrt(P(1,1)))*cos(theta) + a;
      y = (1/sqrt(P(2,2)))*sin(theta) + b;
      plot(x,y,c, 'LineWidth', 3);
      hold on

end 