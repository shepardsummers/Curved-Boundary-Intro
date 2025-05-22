function f_new = GZS_scheme(x1,y1,x2,y2,R,x_circ,y_circ,U_f, U_ff,w,Rho,Ksi,c_s,f,f_eq,ff,ff_eq,Tau)

    % Can be sped up by precomputing w*Rho

    c_w = find_the_wall_point(x1,y1,x2,y2, R, x_circ, y_circ);
    delta = sqrt((c_w(1) - x2)^2 + (c_w(2) - y2)^2)/sqrt((x1-x2)^2 + (y1-y2)^2);

    U_b1 = (delta-1)*U_f/delta;

    if delta >= 0.75
        U_b = U_b1;
        f_a_eq = w*Rho*(1 + U_b'*Ksi/c_s^2 + (U_b'*Ksi)^2/(2*c_s^4) - dot(U_b,U_b)/(2*c_s^2));
        f_a_neq = f - f_eq;
    elseif delta < 0.75
        U_b2 = (delta-1)*U_ff/(1+delta);
        U_b = delta*U_b1 + (1-delta)*U_b2;
        f_a_eq = w*Rho*(1 + U_b'*Ksi/c_s^2 + (U_b'*Ksi)^2/(2*c_s^4) - dot(U_b,U_b)/(2*c_s^2));
        f_a_neq = f - f_eq + (1-delta)*(ff - ff_eq);
    end

    f_new = f_a_eq + (1 - 1/Tau)*f_a_neq;
end

