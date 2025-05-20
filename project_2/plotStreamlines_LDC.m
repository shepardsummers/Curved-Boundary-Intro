for i=1:N_x
    x(i)=(i-1)/(N_x-1);
end
for j=1:N_y
    y(j)=(j-1)/(N_y-1)*((N_y-1)/(N_x-1));
end
%% Generate random seeds
numsamples=floor(0.8*N_x);
seedX=zeros(1,numsamples);
seedY=zeros(1,numsamples);
for i=1:numsamples
    seedX(i)=randi(N_x)/(N_x-1);
end
for i=1:numsamples
    seedY(i)=randi(N_y)/(N_y-1)*((N_y-1)/(N_x-1));
end
%% Generate Fibonacci seeds
Num=55;
[seedX_fibo, seedY_fibo] = fibo(Num, x(1), x(end), y(1), y(end));

%% Combine seeds and plot
seedX=[seedX,seedX_fibo];
seedY=[seedY,seedY_fibo];
plotStreamlines(x, y, flipud(u), flipud(v), seedX, seedY)


function plotStreamlines(x, y, u, v, seedX, seedY)
    % PLOTSTREAMLINES Visualizes streamlines of a flow field.
    % 
    % Inputs:
    %   x     - X-coordinate grid (matrix)
    %   y     - Y-coordinate grid (matrix)
    %   u     - X-component of velocity field (matrix)
    %   v     - Y-component of velocity field (matrix)
    %   seedX - X-coordinates of seed points (vector)
    %   seedY - Y-coordinates of seed points (vector)
    %
    % Example:
    %   [X, Y] = meshgrid(0:0.1:2, 0:0.1:2);
    %   U = -Y; 
    %   V = X;
    %   seedX = [0.5, 1, 1.5]; 
    %   seedY = [0.5, 1, 1.5];
    %   plotStreamlines(X, Y, U, V, seedX, seedY);

    figure;
    hold on;
    % quiver(x, y, u, v, 'r'); % Overlay velocity vectors
    lineobj = streamline(x, y, u, v, seedX, seedY);
    xlabel('X');
    ylabel('Y');
    title('Flow Field Streamlines');
    axis equal;
    L=length(lineobj);
    for i=1:L
        lineobj(i).LineWidth=1.3;
        lineobj(i).Color='black';
    end
    % grid on;
    hold off;
end

function [seedX_fibo, seedY_fibo] = fibo(n, x1, x2, y1, y2)

    % FIBO generates X and Y seeds with a concentration near the four
    % corners with fibonacci series
    % 
    % Inputs:
    %   n     - The number of entries in the Fibonacci series (scalar)
    %   x1    - Lower limit in x (scalar)
    %   x2    - Upper limit in x (scalar)
    %   y1    - Lower limit in y (scalar)
    %   y2    - Upper limit in y (scalar)
    %   seedX_fibo - X-coordinates of seed points (vector)
    %   seedY_fibo - Y-coordinates of seed points (vector)
    %
    s=fibonacci(1:n);

    % x direction
    L=x2-x1;
    L1=L*0.25; % Beginning section -- Fibonacci region
    L2=L*0.5; % Middle section -- Even spaced region
    L3=L*0.25; % Ending section -- Reverse Fibonacci region
    dx1=L1/sum(s)*s;
    dx3=L3/sum(s)*fliplr(s);
    total_num_dx=length(dx1)+floor(L2/dx1(end))+length(dx3);
    total_num_points=total_num_dx+1;
    dx2=L2/floor(L2/dx1(end));
    dx=[dx1,ones(1,floor(L2/dx1(end)))*dx2,dx3];
    % Check
    if single(sum(dx))~=L
        error('dx does not add up');
    end
    seedX_fibo=zeros(1,total_num_points);
    for i=2:total_num_points
        seedX_fibo(i)=seedX_fibo(i-1)+dx(i-1);
    end

    % y direction
    H=y2-y1;
    H1=H*0.25; % Beginning section -- Fibonacci region
    H2=H*0.5; % Middle section -- Even spaced region
    H3=H*0.25; % Ending section -- Reverse Fibonacci region
    dy1=H1/sum(s)*s;
    dy3=H3/sum(s)*fliplr(s);
    total_num_dy=length(dy1)+floor(H2/dy1(end))+length(dy3);
    total_num_points=total_num_dy+1;
    dy2=H2/floor(H2/dy1(end));
    dy=[dy1,ones(1,floor(H2/dy1(end)))*dy2,dy3];
    % Check
    if single(sum(dy))~=H
        error('dx does not add up');
    end
    seedY_fibo=zeros(1,total_num_points);
    for i=2:total_num_points
        seedY_fibo(i)=seedY_fibo(i-1)+dy(i-1);
    end
end