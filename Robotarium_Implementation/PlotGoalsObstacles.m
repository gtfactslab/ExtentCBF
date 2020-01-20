function PlotGoalsObstacles()

%     plot_ellipse(P1, c1, c2, 'b');
% %     dim = [.62 .31 .3 .3];
% %     str = 'Safe Set';
% %     %annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','none');
% %     text(0.5, 0.6, str, 'Color','black','FontSize',14, 'FontWeight', 'bold');
%     hold on   

    p = 4;
    sigma = [1, 0.4, pi/20];
    theta_k = pi/2;
    k = theta_k/(2*sigma(1));
    c = 1;
    theta0 = sign(k)*pi/2;
    [x,y] = meshgrid(-1.6:0.001:1.6, -1:0.001:1);

    x_new = k.*x;
    y_new = k.*y + 1;

    R = ((x_new).^2 + (y_new).^2).^(1/2);
    theta = atan2(y_new, x_new);

    alpha = abs(R - c)/sigma(2);
    beta = abs(theta - theta0)/sigma(1);

    Hg = abs(k) - (alpha.^(p) + beta.^(p)).^(1/p);

    contour(x, y, Hg,[0, 0], 'LineWidth', 4, 'ShowText', 'off');
    hold on
    
end