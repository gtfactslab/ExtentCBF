function h = plot_squircle(x_curr, shape)  
    % Sample the squircle
    th = 0:pi/20:2*pi;
    
    % Parameters
    a = shape(1, 1);
    b = shape(2, 2);
    t = x_curr(3);
    
    rot = [cos(t), sin(t); -sin(t), cos(t)]*[0 1; -1 0];
    holder = [];    
    
    for i = 1:length(th)
        x_unit = sqrt(abs(cos(th(i))))*a*sign(cos(th(i)));
        y_unit = sqrt(abs(sin(th(i))))*b*sign(sin(th(i)));
        holder = [holder; [x_unit, y_unit]*rot + [x_curr(1),x_curr(2)]]; 
    end
    
    h = fill(holder(:, 1), holder(:, 2),[.447,.643,.831],'EdgeColor','none','FaceAlpha',0.7);
end