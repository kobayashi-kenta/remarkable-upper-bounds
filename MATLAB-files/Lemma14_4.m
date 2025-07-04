function Lemma14_4
    pair = [];

    a = sym('a'); at = sym('at'); b = sym('b'); bt = sym('bt');
    d1 = 1 - a; d2 = 1 - 2*a;
    phi = b^2/(a^2 + b^2);
    psi = b^2*(1 + 2*a)/(1 + 4*b^2)*phi;
    psi_a = diff(psi, a);
    psi_b = diff(psi, b);
    psi_aa = diff(psi_a, a);
    psi_bb = diff(psi_b, b);

    f = L1ab(a, b);
    g = (1 + a^2 + (1 - a)^2 + 2*b^2)/28 ...
        - b^4/16/(a^2 + b^2)/((1 - a)^2 + b^2);
    h = (1 - a + a^2 + b^2)/14 ...
        - b^4*(1 + 2*a)/16/(a^2 + b^2)/(1 + 4*b^2) ...
        - b^4*(3 - 2*a)/16/((1 - a)^2 + b^2)/(1 + 4*b^2);
    j = (1 - a + a^2 + b^2)/14 - (psi + subs(psi, a, 1 - a))/16;
    pair = [pair [f; g] [f; h] [f; j]];

    f = (1 - a + a^2 + b^2)/14 - 1/sym(16)*( ...
            b^2*(1 + 2*a)/(1 + 4*b^2) ...
            + b^2*(1 + 2*(1 - a))/(1 + 4*b^2) ...
        );
    g = (1 - a + a^2 + b^2)/14 - b^2/4/(1 + 4*b^2);
    pair = [pair [f; g]];

    f = diff(L1ab(a, b), a);
    g = - (1 - 2*a)/14 - (psi_a - subs(psi_a, a, 1 - a))/16;
    pair = [pair [f; g]];

    f = - (1 - 2*a)/14 - 1/sym(16)*( ...
            2*b^2/(1 + 4*b^2) ...
            + 2*b*(1 + 2*(1 - a))/3/(1 + 4*b^2) ...
        );
    g = - (1 - 2*a)/14 - b*(3 - 2*a + 3*b)/24/(1 + 4*b^2);
    pair = [pair [f; g]];

    f = - (1 - 2*a)/14 - 1/sym(16)*( ...
            -2*b*(1 + 2*a)/3/(1 + 4*b^2) - 2*b^2/(1 + 4*b^2) ...
        );
    g = - (1 - 2*a)/14 + b*(1 + 2*a + 3*b)/24/(1 + 4*b^2);
    pair = [pair [f; g]];

    f = subs(diff(L1ab(a, b), a, 2), a, at);
    g = 1/sym(7) ...
        - (subs(psi_aa, a, at) + subs(psi_aa, a, 1 - at))/16;
    pair = [pair [f; g]];

    f = 1/sym(7) - 1/sym(16)*( ...
            - 2*(3 + 6*at + 4*b)/3/(1 + 4*b^2) ...
            - 2*(3 + 6*(1 - at) + 4*b)/3/(1 + 4*b^2) ...
        );
    g = 1/sym(7) + (3 + 2*b)/6/(1 + 4*b^2);
    pair = [pair [f; g]];

    f = diff(L1ab(a, b), b);
    g = b/7 - (psi_b + subs(psi_b, a, 1 - a))/16;
    pair = [pair [f; g]];

    f = b/7 - 1/sym(16)*( ...
            (1 + 2*a)*(2 + 16*b - 9*b^2)/10/(1 + 4*b^2) ...
            + (1 + 2*(1 - a))*(2 + 16*b - 9*b^2) ...
                /10/(1 + 4*b^2) ...
        );
    g = b/7 - (2 + 16*b - 9*b^2)/40/(1 + 4*b^2);
    pair = [pair [f; g]];

    f = subs(diff(L1ab(a, b), b, 2), b, bt);
    g = 1/sym(7) ...
        - subs(psi_bb + subs(psi_bb, a, 1 - a), b, bt)/16;
    pair = [pair [f; g]];

    f = 1/sym(7) - 1/sym(16)*( ...
            - 21*(1 + 2*a)/11/(1 + 4*bt^2) ...
            - 21*(1 + 2*(1 - a))/11/(1 + 4*bt^2) ...
        );
    g = 1/sym(7) + 21/sym(44)/(1 + 4*bt^2);
    pair = [pair [f; g]];

    f = 168*(1 + 4*b^2)*( ...
            2*((1 - a + a^2 + b^2)/14 - b^2/4/(1 + 4*b^2)) ...
            + b*( ...
                - (1 - 2*a)/14 ...
                - b*(3 - 2*a + 3*b)/24/(1 + 4*b^2) ...
            ) ...
        );
    g = 6*(d2^2 + 4*a*b)*(1 + 4*b^2) ...
        + b^2*(3 + 14*a + 3*b + 46*(1 - b)^2) ...
        + 2*(3 - b - 5*b^2)^2;
    pair = [pair [f; g]];

    f = 168*(1 + 4*b^2)*( ...
            2*((1 - a + a^2 + b^2)/14 - b^2/4/(1 + 4*b^2)) ...
            - b*( ...
                - (1 - 2*a)/14 ...
                + b*(1 + 2*a + 3*b)/24/(1 + 4*b^2) ...
            ) ...
        );
    g = 6*(d2^2 + 4*d1*b)*(1 + 4*b^2) ...
        + b^2*(3 + 14*d1 + 3*b + 46*(1 - b)^2) ...
        + 2*(3 - b - 5*b^2)^2;
    pair = [pair [f; g]];

    f = 84*(1 + 4*b^2)*( ...
            5*((1 - a + a^2 + b^2)/14 - b^2/4/(1 + 4*b^2)) ...
            - b^2*(1/sym(7) + (3 + 2*b)/6/(1 + 4*b^2)) ...
        );
    g = 2*a*d1 + 2*d2^2*(4 + 15*b^2) ...
        + 22*(1 + 2*b)*(1 - b)^2 + 9*b^2*(1 + 2*(1 - 2*b)^2);
    pair = [pair [f; g]];

    f = 280*(1 + 4*b^2)*( ...
            2*((1 - a + a^2 + b^2)/14 - b^2/4/(1 + 4*b^2)) ...
            + b*(b/7 - (2 + 16*b - 9*b^2)/40/(1 + 4*b^2)) ...
        );
    g = 64*a*d1*b*(1 - 2*b)^2 ...
        + 2*(2*a*d1*(30 + b^2) + d2^2*(20 + 13*b))*(1 - b) ...
        + b^2*(11 + 3*d2^2 + 320*b^2 + 63*d2^2*b);
    pair = [pair [f; g]];

    f = 28*(1 + 4*b^2)*( ...
            2*((1 - a + a^2 + b^2)/14 - b^2/4/(1 + 4*b^2)) ...
            - b*b/7 ...
        );
    g = 1 + 2*(1 - b^2) + d2^2*(1 + 4*b^2);
    pair = [pair [f; g]];

    f = 14*(1 + 4*b^2)*( ...
            4*((1 - a + a^2 + b^2)/14 - b^2/4/(1 + 4*b^2)) ...
            - b^2*(1/sym(7) + 1/sym(2)/(1 + 4*b^2)) ...
        );
    g = 1 + b^2 + d2^2*(1 + 4*b^2) + 2*(1 - 2*b^2)^2;
    pair = [pair [f; g]];

    for p = pair
        if simplifyFraction(p(1) - p(2)) ~= 0
            error('There is something wrong.');
        end
    end
    disp('It is all right.');
end
