function h = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit,'Color',[0.72, 0.506, 0.125], ...
                      'LineWidth', 2, ...
                       'HandleVisibility', 'off');
hold off