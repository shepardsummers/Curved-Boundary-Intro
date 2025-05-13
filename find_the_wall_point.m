function [x, y] = find_the_wall_point(x1, y1, x2, y2, R, center_x, center_y)
    
    m = (y2-y1)/(x2-x1);

    b = (y1) - m*(x1);

    A = m^2+1;
    B = 2*(m*b - center_x - m*center_y);
    C = center_x^2 + center_y^2 +b^2 - 2*center_y*b - R^2;

    xp1 = (-B + sqrt(B^2 - 4*A*C))/(2*A);

    xp2 = (-B - sqrt(B^2 - 4*A*C))/(2*A);
    
    if (xp1 >= x1 && xp1 <= x2) || (xp1 <= x1 && xp1 >= x2)
        x = xp1;
    else
        x = xp2;
    end
    
    y = m*x+b;
end

