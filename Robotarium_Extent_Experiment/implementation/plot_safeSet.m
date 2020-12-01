function plot_safeSet(P, c)
      theta = 0:0.001:2*pi;
      x = (1/sqrt(P(1,1)))*cos(theta);
      y = (1/sqrt(P(2,2)))*sin(theta);
      plot(x,y,c, 'LineWidth', 3);
      hold on
end 