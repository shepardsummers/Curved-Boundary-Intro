function s = test_circle(x, y, R, center_x, center_y)
%%% Athis function tests if the point [x;y] is within the circle whose
%%% center's coordinates are [center_x; center_y] and the radius of the
%%% circle is R.

distance = sqrt((x-center_x)^2 + (y-center_y)^2);

if single(distance) <= single(R)
    s = 1;
else
    s = 0;
end