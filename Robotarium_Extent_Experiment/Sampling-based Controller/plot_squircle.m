function h = plot_squircle(x_curr)    
    global shape
    
    th = 0:pi/20:2*pi;
    a = 0.1*shape(1, 1);
    b = 0.1*shape(2, 2);
    t = x_curr(3);
    rot = [cos(t), sin(t); -sin(t), cos(t)]*[0 1; -1 0];
    holder = [];
    x_unit = [];
    y_unit = [];
    for i = 1:length(th)
        x_unit(i) = sqrt(abs(cos(th(i))))*a*sign(cos(th(i)));
        y_unit(i) = sqrt(abs(sin(th(i))))*b*sign(sin(th(i)));
    end

    for i = 1:size(th, 2)
       holder = [holder; [x_unit(i), y_unit(i)]*rot + [x_curr(1),x_curr(2)]]; 
    end 

    h = fill(holder(:, 1), holder(:, 2),[.447,.643,.831],'EdgeColor','none','FaceAlpha',0.3);
end