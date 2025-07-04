function Lemma14_9
    a = sym('a'); b = sym('b');
    d1 = 1 - a; d2 = 1 - 2*a;

    S = b/2;
    c1 = 2*((1-a)^2+b^2)/b; c2 = 2*(a^2+b^2)/b; c3 = 2/b;
    u1 = sym('u1'); u2 = sym('u2'); u3 = sym('u3');
    w1 = sym('w1'); w2 = sym('w2'); w3 = sym('w3');
    w4 = sym('w4'); w5 = sym('w5'); w6 = sym('w6');
    w7 = sym('w7'); w8 = sym('w8'); w9 = sym('w9');

    G1 =  F_beta_1(S/4, c1, c2, c3,  0, u3, u2, w7, w5, w3) ...
        + F_beta_1(S/4, c1, c2, c3, u3,  0, u1, w1, w8, w6) ...
        + F_beta_1(S/4, c1, c2, c3, u2, u1,  0, w4, w2, w9) ...
        + F_beta_1(S/4, c1, c2, c3, u1, u2, u3,-w7,-w8,-w9);

    G2 =  F_beta_2(S/4, c1, c2, c3,  0, u3, u2, w7, w5, w3) ...
        + F_beta_2(S/4, c1, c2, c3, u3,  0, u1, w1, w8, w6) ...
        + F_beta_2(S/4, c1, c2, c3, u2, u1,  0, w4, w2, w9) ...
        + F_beta_2(S/4, c1, c2, c3, u1, u2, u3,-w7,-w8,-w9);

    D = c1*c2*c3/16 - 41*(c1 + c2 + c3)/1080 ...
            - (1/c1 + 1/c2 + 1/c3)/5;
    H = 48*c1^2*c2^2*c3^2*(S*D*G2 - G1);

    E = 48*D - c1 - c2 - c3;
    q1 = 2*E + 2*c1 + c2 + c3;
    q2 = 2*E + 2*c2 + c3 + c1;
    q3 = 2*E + 2*c3 + c1 + c2;
    r1 = c1*(12*D*(E + c1) + 1)/(2*E + c1 + 2*c2 + 2*c3)^2;
    r2 = c2*(12*D*(E + c2) + 1)/(2*E + c2 + 2*c3 + 2*c1)^2;
    r3 = c3*(12*D*(E + c3) + 1)/(2*E + c3 + 2*c1 + 2*c2)^2;
    v1 = c2*c3*(w1 + w4);
    v2 = c3*c1*(w2 + w5);
    v3 = c1*c2*(w3 + w6);
    v4 = w1 + w2 - w3 + w4 - w5 + w6 + 2*(w8 + w9);
    v5 = w2 + w3 - w1 + w5 - w6 + w4 + 2*(w9 + w7);
    v6 = w3 + w1 - w2 + w6 - w4 + w5 + 2*(w7 + w8);
    v7 = c1*c2*c3*(4*q1*(w4 - w8 + w9) - q1*(v5 - v6) ...
                - (c2 - c3)*v4);
    v8 = c1*c2*c3*(4*q2*(w5 - w9 + w7) - q2*(v6 - v4) ...
                - (c3 - c1)*v5);
    v9 = c1*c2*c3*(4*q3*(w6 - w7 + w8) - q3*(v4 - v5) ...
                - (c1 - c2)*v6);

    H2 = 16*( ...
          1/c2/c3/q1*(c1*c2*c3*(q1*(u1 - u2 - u3) ...
                - 6*D*(c1 - c2 - c3)*v4) - 96*D*v1)^2 ...
        + 1/c3/c1/q2*(c1*c2*c3*(q2*(u2 - u3 - u1) ...
                - 6*D*(c2 - c3 - c1)*v5) - 96*D*v2)^2 ...
        + 1/c1/c2/q3*(c1*c2*c3*(q3*(u3 - u1 - u2) ...
                - 6*D*(c3 - c1 - c2)*v6) - 96*D*v3)^2 ...
    ) ...
    + ( ...
          1/c2/c3/q1*(q1*(2*c1*v1 - c2*v2 + c3*v3) ...
                + (c2 - c3)*(96*D - c1)*v1 - v7)^2 ...
        + 1/c3/c1/q2*(q2*(2*c2*v2 - c3*v3 + c1*v1) ...
                + (c3 - c1)*(96*D - c2)*v2 - v8)^2 ...
        + 1/c1/c2/q3*(q3*(2*c3*v3 - c1*v1 + c2*v2) ...
                + (c1 - c2)*(96*D - c3)*v3 - v9)^2 ...
    ) ...
    + 4*( ...
          c1/q1/r1*(v1 + r1*c2*c3*(96*D - c1)*v4)^2 ...
        + c2/q2/r2*(v2 + r2*c3*c1*(96*D - c2)*v5)^2 ...
        + c3/q3/r3*(v3 + r3*c1*c2*(96*D - c3)*v6)^2 ...
    ) ...
    + 8*( ...
          (q1*c1 + 3*D*(c2 + c3)*(c2 - c3)^2) ...
                /q1/(12*D*(E + c1) + 1)*v1^2 ...
        + (q2*c2 + 3*D*(c3 + c1)*(c3 - c1)^2) ...
                /q2/(12*D*(E + c2) + 1)*v2^2 ...
        + (q3*c3 + 3*D*(c1 + c2)*(c1 - c2)^2) ...
                /q3/(12*D*(E + c3) + 1)*v3^2 ...
    ) ...
    + 2/E*( ...
          ((c2 - c3)^2 + 32)/(E + c1) ...
                *(c1 + E/(12*D*(E + c1) + 1))*v1^2 ...
        + ((c3 - c1)^2 + 32)/(E + c2) ...
                *(c2 + E/(12*D*(E + c2) + 1))*v2^2 ...
        + ((c1 - c2)^2 + 32)/(E + c3) ...
                *(c3 + E/(12*D*(E + c3) + 1))*v3^2 ...
    ) ...
    + 2/E*( ...
          (8*E*(3*D*c1 - 1)*c1 - (c2 - c3)^2 - 32)*v1^2 ...
        + (8*E*(3*D*c2 - 1)*c2 - (c3 - c1)^2 - 32)*v2^2 ...
        + (8*E*(3*D*c3 - 1)*c3 - (c1 - c2)^2 - 32)*v3^2 ...
        + 48*D*E*( ...
              (c1*c2 - 4)*v1*v2 ...
            + (c2*c3 - 4)*v2*v3 ...
            + (c3*c1 - 4)*v3*v1 ...
        ) ...
    );

    e1 = 45*b^3*E/4;
    e2 = 108*(a*d1 - b^2)^2*(1 + d2^2 + 4*b^2)/c1/c2 ...
        + 270*a^2*d1^2 ...
        + b^2*(47*a*d1 + 88*d2^2 + (1 - 100*b^2)) + 81*b^4;
    if simplifyFraction(e1 - e2) ~= 0
        error('There is something wrong.');
    end

    %%  Verify the equality term by term.  %%
    variables = [u1 u2 u3 w1 w2 w3 w4 w5 w6 w7 w8 w9];
    for p = 1:12
        df1 = diff(H - H2, variables(p));
        for q = p:12
            df2 = simplifyFraction(diff(df1, variables(q)));
            if df2 ~= 0
                error('There is something wrong.');
            end
        end
    end
    disp('It is all right.');
end
