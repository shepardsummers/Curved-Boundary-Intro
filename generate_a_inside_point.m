function [x, y] = generate_a_inside_point(R, center_x, center_y)

    [x,y] = generate_a_point();

    while R < sqrt((x-center_x)^2+(y-center_y)^2)
        [x, y] = generate_a_point();
    end

end

