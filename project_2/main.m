clear all
%% Defination of Parameters
% Domain Related
D = 21; 

N_x = 35*D;
N_y = 9*D;
dx=1;
dy=1;

N_x_circ_center = round(5*D);
N_y_circ_center = round(N_y/2);
center_x = dx * (N_x_circ_center - 1);
center_y = dy * (N_y_circ_center - 1);
R = D/2;

for i = 1:N_x
    x(i)=dx*(i-1);
end
for j = 1:N_y
    y(j)=dy*(j-1);
end
% Domain ID
% Domain=0; Solid
% Domain=1; Fluid
for j = 1:N_y
    for i = 1:N_x
        if test_circle(x(i), y(j), R, center_x, center_y)
            Domain_ID(j,i) = 0;
        else 
            Domain_ID(j,i) = 1;
        end
    end
end
%figure
%contourf(Domain_ID, 30)
%axis equal tight

% Zone ID
% Zone ID=0 --> dead zone
% Zone_ID=1 --> b nodes
% Zone_ID=2 --> f nodes
% Zone_ID=3 --> other nodes in fluid domain
for j = 1:N_y
    for i = 1:N_x
        if Domain_ID(j,i) == 0 % solid domain
            if Domain_ID(j,i+1) == 1 || Domain_ID(j,i-1) == 1 || Domain_ID(j+1,i) == 1 || Domain_ID(j-1,i) == 1 ...
                    || Domain_ID(j+1,i+1) == 1 || Domain_ID(j-1,i+1) == 1 || Domain_ID(j+1,i-1) == 1 || Domain_ID(j-1,i-1) == 1
                Zone_ID(j,i) = 1;
            else
                Zone_ID(j,i) = 0;
            end
        else % fluid domain
            if i == 1 || i == N_x || j == 1 || j == N_y
                Zone_ID(j,i) = 3;
            else
                if Domain_ID(j,i+1) == 0 || Domain_ID(j,i-1) == 0 || Domain_ID(j+1,i) == 0 || Domain_ID(j-1,i) == 0 ...
                    || Domain_ID(j+1,i+1) == 0 || Domain_ID(j-1,i+1) == 0 || Domain_ID(j+1,i-1) == 0 || Domain_ID(j-1,i-1) == 0
                    Zone_ID(j,i) = 2;
                else
                    Zone_ID(j,i) = 3;
                end
            end
            
        end
    end
end
figure
contourf(Zone_ID, 30)
axis equal tight

% LBM Related
Ksi=[0 1 0 -1  0 1  -1  -1  1;...
     0 0 1  0 -1 1   1  -1 -1]; % Lattice Velocities for D2Q9
w=[4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36]; % Weights for D2Q9
c_s=1/sqrt(3); % Speed of Sound for D2Q9
Tau=1.2; % Relaxation Time
Rho_in=2;

Re = 20;
viscos = (Tau-0.5)*c_s^2;
U_in = Re*viscos/(D);

%% Initialization
Rho_ref=2;
Rho=ones(1,N_y,N_x)*Rho_ref;
U=zeros(2,N_y,N_x);

f=zeros(9,N_y,N_x); % PDF for all 9 directions and at all locations
for j=1:N_y
    for i=1:N_x
        f(:,j,i) = eqm_d2q9(squeeze(Rho(1,j,i)),squeeze(U(:,j,i)));
    end
end
f_new=f;
f_eq=f;

Timer=50000;
%% Solving
x_circ = center_x;
y_circ = center_y;
for t=1:Timer
% Streaming/Boundary Conditions
for j=1:N_y
    for i=1:N_x
        if Zone_ID(j,i) == 0 % dead zone
            % Do nothin
        elseif Zone_ID(j,i) == 1 % b node
            % Do nothin
        elseif Zone_ID(j,i) == 2 % f node
            % FH curved boundary scheme goes here
            %% Direction 1
            f_new(1,j,i) = f(1,j,i);

            %% Direction 2
            if Zone_ID(j,i+1) == 1 % b-node
                % xy coords of node inside of circle
                x1=x(i+1);
                y1=y(j);
                % xy coords of node outside of circle
                x2=x(i);
                y2=y(j);
                
                c_w = find_the_wall_point(x1,y1,x2,y2, R, x_circ, y_circ);
                delta = sqrt((c_w(1) - x2)^2 + (c_w(2) - y2)^2)/sqrt((x1-x2)^2 + (y1-y2)^2);
                if delta >= 0.5
                    Chi = (2*delta - 1)/Tau;
                    U_bf = (delta - 1)*squeeze(U(:,j,i))/delta;
                elseif delta < 0.5
                    Chi = (2*delta - 1)/(Tau - 1);
                    U_bf = squeeze(U(:,j,i));
                end
                U_f = squeeze(U(:,j,i))';
                
                f_star = w(2)*Rho(1,j,i)*(1 + U_bf'*Ksi(:,2)/c_s^2 + (U_f*Ksi(:,2))^2/(2*c_s^4) - dot(U_f,U_f)/(2*c_s^2));
                
                f_new(4,j,i) = (1 - Chi)*f(2,j,i) + Chi*f_star; % Weird "bounceback"
            else
                f_new(4,j,i) = f(4,j,i+1);
            end

            %% Direction 3
            if Zone_ID(j-1,i) == 1 % b-node
                % xy coords of node inside of circle
                x1=x(i);
                y1=y(j-1);
                % xy coords of node outside of circle
                x2=x(i);
                y2=y(j);
                
                c_w = find_the_wall_point(x1,y1,x2,y2, R, x_circ, y_circ);
                delta = sqrt((c_w(1) - x2)^2 + (c_w(2) - y2)^2)/sqrt((x1-x2)^2 + (y1-y2)^2);
                if delta >= 0.5
                    Chi = (2*delta - 1)/Tau;
                    U_bf = (delta - 1)*squeeze(U(:,j,i))/delta;
                elseif delta < 0.5
                    Chi = (2*delta - 1)/(Tau - 1);
                    U_bf = squeeze(U(:,j,i));
                end
                U_f = squeeze(U(:,j,i))';
                
                f_star = w(3)*Rho(1,j,i)*(1 + U_bf'*Ksi(:,3)/c_s^2 + (U_f*Ksi(:,3))^2/(2*c_s^4) - dot(U_f,U_f)/(2*c_s^2));
                
                f_new(5,j,i) = (1 - Chi)*f(3,j,i) + Chi*f_star; % Weird "bounceback"
            else
                f_new(5,j,i) = f(5,j-1,i);
            end
            %% Direction 4
            if Zone_ID(j,i-1) == 1 % b-node
                % xy coords of node inside of circle
                x1=x(i-1);
                y1=y(j);
                % xy coords of node outside of circle
                x2=x(i);
                y2=y(j);
                
                c_w = find_the_wall_point(x1,y1,x2,y2, R, x_circ, y_circ);
                delta = sqrt((c_w(1) - x2)^2 + (c_w(2) - y2)^2)/sqrt((x1-x2)^2 + (y1-y2)^2);
                if delta >= 0.5
                    Chi = (2*delta - 1)/Tau;
                    U_bf = (delta - 1)*squeeze(U(:,j,i))/delta;
                elseif delta < 0.5
                    Chi = (2*delta - 1)/(Tau - 1);
                    U_bf = squeeze(U(:,j,i));
                end
                U_f = squeeze(U(:,j,i))';
                
                f_star = w(4)*Rho(1,j,i)*(1 + U_bf'*Ksi(:,4)/c_s^2 + (U_f*Ksi(:,4))^2/(2*c_s^4) - dot(U_f,U_f)/(2*c_s^2));
                
                f_new(2,j,i) = (1 - Chi)*f(4,j,i) + Chi*f_star; % Weird "bounceback"
            else
                f_new(2,j,i) = f(2,j,i-1);
            end
            %% Direction 5
            if Zone_ID(j+1,i) == 1 % b-node
                % xy coords of node inside of circle
                x1=x(i);
                y1=y(j+1);
                % xy coords of node outside of circle
                x2=x(i);
                y2=y(j);
                
                c_w = find_the_wall_point(x1,y1,x2,y2, R, x_circ, y_circ);
                delta = sqrt((c_w(1) - x2)^2 + (c_w(2) - y2)^2)/sqrt((x1-x2)^2 + (y1-y2)^2);
                if delta >= 0.5
                    Chi = (2*delta - 1)/Tau;
                    U_bf = (delta - 1)*squeeze(U(:,j,i))/delta;
                elseif delta < 0.5
                    Chi = (2*delta - 1)/(Tau - 1);
                    U_bf = squeeze(U(:,j,i));
                end
                U_f = squeeze(U(:,j,i))';
                
                f_star = w(5)*Rho(1,j,i)*(1 + U_bf'*Ksi(:,5)/c_s^2 + (U_f*Ksi(:,5))^2/(2*c_s^4) - dot(U_f,U_f)/(2*c_s^2));
                
                f_new(3,j,i) = (1 - Chi)*f(5,j,i) + Chi*f_star; % Weird "bounceback"
            else
                f_new(3,j,i) = f(3,j+1,i);
            end
            %% Direction 6
            if Zone_ID(j-1,i+1) == 1 % b-node
                % xy coords of node inside of circle
                x1=x(i+1);
                y1=y(j-1);
                % xy coords of node outside of circle
                x2=x(i);
                y2=y(j);
                
                c_w = find_the_wall_point(x1,y1,x2,y2, R, x_circ, y_circ);
                delta = sqrt((c_w(1) - x2)^2 + (c_w(2) - y2)^2)/sqrt((x1-x2)^2 + (y1-y2)^2);
                if delta >= 0.5
                    Chi = (2*delta - 1)/Tau;
                    U_bf = (delta - 1)*squeeze(U(:,j,i))/delta;
                elseif delta < 0.5
                    Chi = (2*delta - 1)/(Tau - 1);
                    U_bf = squeeze(U(:,j,i));
                end
                U_f = squeeze(U(:,j,i))';
                
                f_star = w(6)*Rho(1,j,i)*(1 + U_bf'*Ksi(:,6)/c_s^2 + (U_f*Ksi(:,6))^2/(2*c_s^4) - dot(U_f,U_f)/(2*c_s^2));
                
                f_new(8,j,i) = (1 - Chi)*f(6,j,i) + Chi*f_star; % Weird "bounceback"
            else
                f_new(8,j,i) = f(8,j-1,i+1);
            end
            %% Direction 7
            if Zone_ID(j-1,i-1) == 1 % b-node
                % xy coords of node inside of circle
                x1=x(i-1);
                y1=y(j-1);
                % xy coords of node outside of circle
                x2=x(i);
                y2=y(j);
                
                c_w = find_the_wall_point(x1,y1,x2,y2, R, x_circ, y_circ);
                delta = sqrt((c_w(1) - x2)^2 + (c_w(2) - y2)^2)/sqrt((x1-x2)^2 + (y1-y2)^2);
                if delta >= 0.5
                    Chi = (2*delta - 1)/Tau;
                    U_bf = (delta - 1)*squeeze(U(:,j,i))/delta;
                elseif delta < 0.5
                    Chi = (2*delta - 1)/(Tau - 1);
                    U_bf = squeeze(U(:,j,i));
                end
                U_f = squeeze(U(:,j,i))';
                
                f_star = w(7)*Rho(1,j,i)*(1 + U_bf'*Ksi(:,7)/c_s^2 + (U_f*Ksi(:,7))^2/(2*c_s^4) - dot(U_f,U_f)/(2*c_s^2));
                
                f_new(9,j,i) = (1 - Chi)*f(7,j,i) + Chi*f_star; % Weird "bounceback"
            else
                f_new(9,j,i) = f(9,j-1,i-1);
            end
            %% Direction 8
            if Zone_ID(j+1,i-1) == 1 % b-node
                % xy coords of node inside of circle
                x1=x(i-1);
                y1=y(j+1);
                % xy coords of node outside of circle
                x2=x(i);
                y2=y(j);
                
                c_w = find_the_wall_point(x1,y1,x2,y2, R, x_circ, y_circ);
                delta = sqrt((c_w(1) - x2)^2 + (c_w(2) - y2)^2)/sqrt((x1-x2)^2 + (y1-y2)^2);
                if delta >= 0.5
                    Chi = (2*delta - 1)/Tau;
                    U_bf = (delta - 1)*squeeze(U(:,j,i))/delta;
                elseif delta < 0.5
                    Chi = (2*delta - 1)/(Tau - 1);
                    U_bf = squeeze(U(:,j,i));
                end
                U_f = squeeze(U(:,j,i))';
                
                f_star = w(8)*Rho(1,j,i)*(1 + U_bf'*Ksi(:,8)/c_s^2 + (U_f*Ksi(:,8))^2/(2*c_s^4) - dot(U_f,U_f)/(2*c_s^2));
                
                f_new(6,j,i) = (1 - Chi)*f(8,j,i) + Chi*f_star; % Weird "bounceback"
            else
                f_new(6,j,i) = f(6,j+1,i-1);
            end
            %% Direction 9
            if Zone_ID(j+1,i+1) == 1 % b-node
                % xy coords of node inside of circle
                x1=x(i+1);
                y1=y(j+1);
                % xy coords of node outside of circle
                x2=x(i);
                y2=y(j);
                
                c_w = find_the_wall_point(x1,y1,x2,y2, R, x_circ, y_circ);
                delta = sqrt((c_w(1) - x2)^2 + (c_w(2) - y2)^2)/sqrt((x1-x2)^2 + (y1-y2)^2);
                if delta >= 0.5
                    Chi = (2*delta - 1)/Tau;
                    U_bf = (delta - 1)*squeeze(U(:,j,i))/delta;
                elseif delta < 0.5
                    Chi = (2*delta - 1)/(Tau - 1);
                    U_bf = squeeze(U(:,j,i));
                end
                U_f = squeeze(U(:,j,i))';
                
                f_star = w(9)*Rho(1,j,i)*(1 + U_bf'*Ksi(:,9)/c_s^2 + (U_f*Ksi(:,9))^2/(2*c_s^4) - dot(U_f,U_f)/(2*c_s^2));
                
                f_new(7,j,i) = (1 - Chi)*f(9,j,i) + Chi*f_star; % Weird "bounceback"
            else
                f_new(7,j,i) = f(7,j+1,i+1);
            end
        else
            % All other nodes in fluid domain
            if j == 1 % This is the top boundary nodes
                if i ==1 % Top-Left corner node
                    f_new(1,j,i) = f(1,j,i);
                    f_new(3,j,i) = f(3,j+1,i);
                    f_new(4,j,i) = f(4,j,i+1);
                    f_new(7,j,i) = f(7,j+1,i+1);
    
                    % Unknown
                    f_new(2,j,i) = f_new(4,j,i)+Rho_in*U_in*2/3;
                    f_new(5,j,i) = f_new(5,N_y,i);
                    %f_new(6,j,i) = f_new(6,N_y,i);
                    f_new(8,j,i) = f_new(8,N_y,i);
                    %f_new(9,j,i) = f_new(9,N_y,i);
                    f_new(6,j,i) = f_new(8,j,i)+(f_new(5,j,i)-f_new(3,j,i))/2+Rho_in*U_in/6; % Double Check
                    f_new(9,j,i) = f_new(7,j,i)-(f_new(5,j,i)-f_new(3,j,i))/2+Rho_in*U_in/6; % Double Check
                elseif i == N_x % Top-right corner node
                    f_new(1,j,i) = f(1,j,i);
                    f_new(2,j,i) = f(2,j,i-1);
                    f_new(3,j,i) = f(3,j+1,i);
                    f_new(6,j,i) = f(6,j+1,i-1);
    
                    % Unknown
                    %f_new(4,j,i) = f_new(2,j,i);
                    f_new(5,j,i) = f_new(5,N_y,i);
                    %f_new(7,j,i) = f_new(7,j,i-1);
                    %f_new(8,j,i) = f_new(6,j,i);
                    f_new(9,j,i) = f_new(9,N_y,i);
                    f_new(4,j,i) = f_new(4,j,i-1);
                    f_new(7,j,i) = f_new(7,j,i-1);
                    f_new(8,j,i) = f_new(8,j,i-1);
                else % All other nodes on the top boundary
                    f_new(1,j,i) = f(1,j,i);
                    f_new(2,j,i) = f(2,j,i-1);
                    f_new(3,j,i) = f(3,j+1,i);
                    f_new(4,j,i) = f(4,j,i+1);
                    f_new(6,j,i) = f(6,j+1,i-1);
                    f_new(7,j,i) = f(7,j+1,i+1);
    
                    % Unknown
                    f_new(5,j,i) = f_new(5,N_y,i);
                    f_new(8,j,i) = f_new(8,N_y,i);
                    f_new(9,j,i) = f_new(9,N_y,i);
                end
            elseif j == N_y % This is the bottom boundary nodes
                if i ==1 % Bottom-Left corner node
                    f_new(1,j,i) = f(1,j,i);
                    f_new(4,j,i) = f(4,j,i+1);
                    f_new(5,j,i) = f(5,j-1,i);
                    f_new(8,j,i) = f(8,j-1,i+1);
                    
                    % Unknown
                    f_new(2,j,i) = f_new(4,j,i)+Rho_in*U_in*2/3;
                    f_new(3,j,i) = f_new(3,1,i);
                    %f_new(6,j,i) = f_new(6,1,i);
                    f_new(7,j,i) = f_new(7,1,i);
                    %(9,j,i) = f_new(9,1,i);
                    f_new(6,j,i) = f_new(8,j,i)+(f_new(5,j,i)-f_new(3,j,i))/2+Rho_in*U_in/6; % Double Check
                    f_new(9,j,i) = f_new(7,j,i)-(f_new(5,j,i)-f_new(3,j,i))/2+Rho_in*U_in/6; % Double Check
                elseif i == N_x % Bottom-right corner node
                    f_new(1,j,i) = f(1,j,i);
                    f_new(2,j,i) = f(2,j,i-1);
                    f_new(5,j,i) = f(5,j-1,i);
                    f_new(9,j,i) = f(9,j-1,i-1);
    
                    % Unknown
                    f_new(4,j,i) = f_new(4,j,i-1);
                    f_new(7,j,i) = f_new(7,j,i-1);
                    f_new(8,j,i) = f_new(8,j,i-1);
                    f_new(3,j,i) = f_new(3,1,i);
                    %f_new(4,j,i) = f_new(2,j,i);
                    f_new(6,j,i) = f_new(6,1,i);
                    %f_new(7,j,i) = f_new(9,j,i);
                    %f_new(8,j,i) = f_new(6,j,i);
                else % All other nodes on the bottom boundary
                    f_new(1,j,i) = f(1,j,i);
                    f_new(2,j,i) = f(2,j,i-1);
                    f_new(4,j,i) = f(4,j,i+1);
                    f_new(5,j,i) = f(5,j-1,i);
                    f_new(8,j,i) = f(8,j-1,i+1);
                    f_new(9,j,i) = f(9,j-1,i-1);
    
                    % Unknown
                    f_new(3,j,i) = f_new(3,1,i);
                    f_new(6,j,i) = f_new(6,1,i);
                    f_new(7,j,i) = f_new(7,1,i);
                end
            elseif i == 1 % This is the left boundary nodes
                f_new(1,j,i) = f(1,j,i);
                f_new(3,j,i) = f(3,j+1,i);
                f_new(4,j,i) = f(4,j,i+1);
                f_new(5,j,i) = f(5,j-1,i);
                f_new(7,j,i) = f(7,j+1,i+1);
                f_new(8,j,i) = f(8,j-1,i+1);
    
                % Unknown
                %U_in = 1-(f_new(1,j,i)+f_new(3,j,i)+f_new(5,j,i)+2*(f_new(4,j,i)+f_new(7,j,i)+f_new(8,j,i)))/Rho_in;
                f_new(2,j,i) = f_new(4,j,i)+Rho_in*U_in*2/3;
                f_new(6,j,i) = f_new(8,j,i)+(f_new(5,j,i)-f_new(3,j,i))/2+Rho_in*U_in/6; % Double Check
                f_new(9,j,i) = f_new(7,j,i)-(f_new(5,j,i)-f_new(3,j,i))/2+Rho_in*U_in/6; % Double Check
            elseif i == N_x % This is the right boundary nodes
                f_new(1,j,i) = f(1,j,i);
                f_new(2,j,i) = f(2,j,i-1);
                f_new(3,j,i) = f(3,j+1,i);
                f_new(5,j,i) = f(5,j-1,i);
                f_new(6,j,i) = f(6,j+1,i-1);
                f_new(9,j,i) = f(9,j-1,i-1);
    
                % Unknown
                f_new(4,j,i) = f_new(4,j,i-1);
                f_new(7,j,i) = f_new(7,j,i-1);
                f_new(8,j,i) = f_new(8,j,i-1);
            else  % All interior nodes
                f_new(1,j,i) = f(1,j,i);
                f_new(2,j,i) = f(2,j,i-1);
                f_new(3,j,i) = f(3,j+1,i);
                f_new(4,j,i) = f(4,j,i+1);
                f_new(5,j,i) = f(5,j-1,i);
                f_new(6,j,i) = f(6,j+1,i-1);
                f_new(7,j,i) = f(7,j+1,i+1);
                f_new(8,j,i) = f(8,j-1,i+1);
                f_new(9,j,i) = f(9,j-1,i-1);
            end
        end
    end
end

old_rho = sum(Rho, "all");

% Collision
% Rho, U calculation
Rho = sum(f_new, 1);
U = pagemtimes(Ksi,f_new)./Rho;

% f_eq calculation
% for j=1:N_y
%     for i=1:N_x
%         f_eq(:,j,i) = eqm_d2q9(squeeze(Rho(1,j,i)),squeeze(U(:,j,i)));
%     end
% end
f_eq = pagemtimes(w',Rho).*(1 + 3*pagemtimes(Ksi',U) + 9/2*(pagemtimes(Ksi', U).^2) - 3/2*sum(U.*U,1));
% BGK Collision & Update
f=f_new-(f_new-f_eq)/Tau;
err = sum(Rho,"all")-old_rho;

fprintf("Itt: %i   |   Err: %e\n", t, err)
end
%% Post-Processing/Visualizaation
y_benchmark=0:0.01:1;
for i=1:length(y_benchmark)
    u_benchmark(i)=-4*(y_benchmark(i)^2-y_benchmark(i));
end
figure
plot(y_benchmark,u_benchmark,"red")

% Sampling u_x velocity from the outlet
for j=1:N_y
    u_sim(j) = U(1,j,N_x);
end
hold on
plot((0:1:N_y-1)/(N_y-1),u_sim/max(u_sim),"blue")

figure
quiver(flipud(squeeze(U(1,:,:))),flipud(squeeze(U(2,:,:))),10)
axis equal tight

figure
contourf(flipud(squeeze(Rho)),30)
axis equal tight

vel_mag = zeros(N_y, N_x);
for j = 1:N_y
    for i = 1:N_x
        vel_mag(j,i) = sqrt(U(1,j,i)^2 + U(2,j,i)^2);
    end
end
figure
contourf(flipud(vel_mag),30)
axis equal tight

figure
siz = 6;
u = flip(squeeze(U(1, :, :)));
v = squeeze(U(2, :, :));
[startX, startY] = meshgrid(1:siz:N_x, 1:4:N_y);
verts = stream2(1:N_x,1:N_y,u,v,startX,startY);
streamline(verts)
axis equal tight