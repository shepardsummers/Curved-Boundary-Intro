clear; clc;

num = 500000;

%[x, y] = generate_a_point();

R = 0.4;
center_x = 0.5;
center_y = 0.5;

%s = test_circle(x, y, R, center_x, center_y);

[x1, y1] = generate_a_inside_point(R, center_x, center_y);
[x2, y2] = generate_a_outside_point(R, center_y, center_y);

[x, y] = find_the_wall_point(x1, y1, x2, y2, R, center_x, center_y);

%%
figure
hold on
th = 0:pi/720:2*pi;
xunit = R*cos(th) + center_x;
yunit = R*sin(th) + center_y;
h = plot(xunit, yunit, 'Color', "black");
hold on
plot(x1,y1, '+', 'Color', 'blue', LineWidth=1)
hold on
plot(x2,y2, '+', 'Color', 'blue', LineWidth=1)
hold on
plot(x,y, '+', 'Color', 'red', LineWidth=1)
hold off
axis equal tight