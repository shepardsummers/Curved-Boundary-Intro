function out = test_circle(x,y,R,center_x,center_y)
    if R == sqrt((x-center_x)^2+(y-center_y)^2)
        out = true;
    else
        out = false;
    end
end

