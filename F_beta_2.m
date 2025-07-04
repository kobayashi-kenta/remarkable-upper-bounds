function z = F_beta_2(s, c1, c2, c3, u1, u2, u3, w1, w2, w3)
    v4 = 2*u1 - u2 - u3 + (4*w1 - (c2 - c3)*(u2 - u3)) / c1;
    v5 = 2*u2 - u3 - u1 + (4*w2 - (c3 - c1)*(u3 - u1)) / c2;
    v6 = 2*u3 - u1 - u2 + (4*w3 - (c1 - c2)*(u1 - u2)) / c3;
    z = ( ...
        (c1*v4 + c2*v5 + c3*v6)^2 - 8*(v4*v5 + v5*v6 + v6*v4) ...
    ) / s / 16;
end
