function fc = plot_squircle(x_curr, Vsym, x1, x2, theta)

    V = subs(Vsym,{x1, x2, theta}, {x_curr(1), x_curr(2), x_curr(3)});

    fc=fcontour(matlabFunction(V),[-1.2,1.2], 'HandleVisibility', 'off');
    fc.LevelList=[0];
    fc.LineColor = 'k';
    fc.LineWidth = 3;
    drawnow
    
end
