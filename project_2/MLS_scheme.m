function f_new = MLS_scheme(x1,y1,x2,y2,R,x_circ,y_circ,U,w,Rho,Ksi,c_s,f,Tau,U_ff)

    c_w = find_the_wall_point(x1,y1,x2,y2, R, x_circ, y_circ);
    delta = sqrt((c_w(1) - x2)^2 + (c_w(2) - y2)^2)/sqrt((x1-x2)^2 + (y1-y2)^2);
    if delta >= 0.5
        Chi = (2*delta - 1)/Tau;
        U_bf = (delta - 1)*U/delta;
    elseif delta < 0.5
        Chi = (2*delta - 1)/(Tau - 2);
        U_bf = U_ff;
    end
    U_f = U';
    
    f_star = w*Rho*(1 + U_bf'*Ksi/c_s^2 + (U_f*Ksi)^2/(2*c_s^4) - dot(U_f,U_f)/(2*c_s^2));
    
    f_new = (1 - Chi)*f + Chi*f_star; % Weird "bounceback"
end

