function h = plot_ellipse(P, a, b, phi)

      theta = 0:0.1:2*pi;
      x = (1/sqrt(P(1,1)))*cos(theta);
      y = (1/sqrt(P(2,2)))*sin(theta);
      
      rot = [cos(phi) -sin(phi); sin(phi) cos(phi)];
      holder = rot*[x; y];
      h = plot(holder(1,:)+a,holder(2,:)+b, 'Color',[0.72, 0.506, 0.125], 'LineWidth', 3);
      hold on

%     a = 0.2;
%     b = 0.1;
%     t = phi;
%     rot = [cos(t), sin(t); -sin(t), cos(t)];
% 
%     holder = [];
%     x_unit = [];
%     y_unit = [];
%     th = 0:pi/20:2*pi;
%     
%     for i = 1:length(th)
%         x_unit(i) = sqrt(abs(cos(th(i))))*a*sign(cos(th(i)));
%         y_unit(i) = sqrt(abs(sin(th(i))))*b*sign(sin(th(i)));
%     end
% 
%     for i = 1:size(th, 2)
%        holder = [holder; [x_unit(i), y_unit(i)]*rot + [x,y]]; 
%     end 
%     h = fill(holder(:, 1), holder(:, 2),[.447,.643,.831],'EdgeColor','none','FaceAlpha',0.3);

end 