function Lemma14_7
    pair = [];

    a = sym('a'); at = sym('at'); b = sym('b'); bt = sym('bt');
    d1 = 1 - a; d2 = 1 - 2*a;
    phi = b^2/(a^2 + b^2);
    phi_a = diff(phi, a);
    phi_b = diff(phi, b);
    phi_aa = diff(phi_a, a);
    phi_bb = diff(phi_b, b);

    f = L4ab(a, b);
    g = (a^2 + b^2)*((1 - a)^2 + b^2)/4/b^2 ...
        - (1 + a^2 + (1 - a)^2 + 2*b^2)/30 ...
        - b^2/20*(1 + 1/(a^2 + b^2) + 1/((1-a)^2 + b^2));
    h = a^2*(1 - a)^2/4/b^2 + (11 - 26*a*(1 - a))/60 ...
        + 2*b^2/15 - (phi + subs(phi, a, 1 - a))/20;
    pair = [pair [f; g] [f; h]];

    f = a^2*(1 - a)^2/4/b^2 + (11 - 26*a*(1 - a))/60 ...
        + 2*b^2/15 - (1 + 1)/20;
    g = a^2*(1 - a)^2/4/b^2 + (5 - 26*a*(1 - a) + 8*b^2)/60;
    pair = [pair [f; g]];

    f = diff(L4ab(a, b), a);
    g = (1 - 2*a)*(a*(1 - a)/2/b^2 - 13/sym(30)) ...
        - (phi_a - subs(phi_a, a, 1 - a))/20;
    pair = [pair [f; g]];

    f = subs(diff(L4ab(a, b), a, 2), a, at);
    g = (1 - 6*at*(1 - at))/2/b^2 + 13/sym(15) ...
        - (subs(phi_aa, a, at) + subs(phi_aa, a, 1 - at))/20;
    pair = [pair [f; g]];

    f = (1 - 6*at*(1 - at))/2/b^2 + 13/sym(15) ...
        + 1/sym(20)*(2/b^2 + 2/b^2);
    g = (7 - 30*a*(1 - a) - 30*(at - a)*(1 - at - a))/10/b^2 ...
        + 13/sym(15);
    pair = [pair [f; g]];

    f = (7 - 30*a*(1 - a) + 30*(b/45))/10/b^2 + 13/sym(15);
    g = (21 - 90*a*(1 - a) + 2*b)/30/b^2 + 13/sym(15);
    pair = [pair [f; g]];

    f = diff(L4ab(a, b), b);
    g = - a^2*(1 - a)^2/2/b^3 + 4*b/15 ...
        - (phi_b + subs(phi_b, a, 1 - a))/20;
    pair = [pair [f; g]];

    f = - a^2*(1 - a)^2/2/b^3 + 4*b/15 ...
        - 1/sym(20)*(1/2/b + 1/2/b);
    g = - a^2*(1 - a)^2/2/b^3 - 1/sym(20)/b + 4*b/15;
    pair = [pair [f; g]];

    f = subs(diff(L4ab(a, b), b, 2), b, bt);
    g = 3*a^2*(1 - a)^2/2/bt^4 + 4/sym(15) ...
        - subs(phi_bb + subs(phi_bb, a, 1 - a), b, bt)/20;
    pair = [pair [f; g]];

    f = 3*a^2*(1 - a)^2/2/bt^4 + 4/sym(15) ...
        + 1/sym(20)*(7/9/bt^2 + 7/9/bt^2);
    g = 3*a^2*(1 - a)^2/2/bt^4 + 7/sym(90)/bt^2 + 4/sym(15);
    pair = [pair [f; g]];

    f = 60*b^2*( ...
            3*(a^2*(1 - a)^2/4/b^2 ...
                + (5 - 26*a*(1 - a))/60 + 2*b^2/15) ...
            + b*((1 - 2*a)*(a*(1-a)/2/b^2 - 13/sym(30)) ...
                - 1/sym(30)/b) ...
        );
    g = 11*a*d1^2*b + 5*(3*a*d1 - b)^2 + a*b*(7*d1 - 5*b)^2 ...
        + b^2*(2*(1 - b) + a*(4 + 2*a + 3*b) + 6*(d1 - 2*b)^2);
    pair = [pair [f; g]];

    f = 60*b^2*( ...
            3*(a^2*(1 - a)^2/4/b^2 ...
                + (5 - 26*a*(1 - a))/60 + 2*b^2/15) ...
            - b*((1 - 2*a)*(a*(1 - a)/2/b^2 - 13/sym(30)) ...
                + 1/sym(30)/b) ...
        );
    g = 11*a^2*d1*b + 5*(3*a*d1 - b)^2 + d1*b*(7*a - 5*b)^2 ...
        + b^2*(2*(1 - b) + d1*(4 + 2*d1 + 3*b) + 6*(a - 2*b)^2);
    pair = [pair [f; g]];

    f = 60*b^2*( ...
            9*(a^2*(1 - a)^2/4/b^2 ...
                + (5 - 26*a*(1 - a))/60 + 2*b^2/15) ...
            - b^2*((21 - 90*a*(1-a) + 2*b)/30/b^2 + 13/sym(15))...
        );
    g = 86*a^2*d1^2 + (7*a*d1 - 4*b^2)^2 ...
        + b^2*(2 + 2*a*d1 + (1 - 2*b)^2);
    pair = [pair [f; g]];

    f = 60*b^2*( ...
            3*(a^2*(1 - a)^2/4/b^2 ...
                + (5 - 26*a*(1 - a))/60 + 2*b^2/15) ...
            + b*(- a^2*(1 - a)^2/2/b^3 - 1/sym(20)/b + 4*b/15) ...
        );
    g = b^2*(12*d2^2 + 25*b^2) + 15*(a*d1 - b^2)^2;
    pair = [pair [f; g]];

    f = 60*b^2*( ...
            3*(a^2*(1 - a)^2/4/b^2 ...
                + (5 - 26*a*(1 - a))/60 + 2*b^2/15) ...
            - b*(- a^2*(1 - a)^2/2/b^3 + 4*b/15) ...
        );
    g = b^2*(3 + 12*d2^2 + 5*b^2) + 3*(5*a*d1 - b^2)^2;
    pair = [pair [f; g]];

    f = 60*b^2*( ...
            9*(a^2*(1 - a)^2/4/b^2 ...
                + (5 - 26*a*(1 - a))/60 + 2*b^2/15) ...
            - b^2*(49*a^2*(1 - a)^2/30/b^4 ...
                + 1/sym(12)/b^2 + 4/sym(15)) ...
        );
    g = b^2*(40*d2^2 + 19*b^2) + 37*(a*d1 - b^2)^2;
    pair = [pair [f; g]];

    for p = pair
        if simplifyFraction(p(1) - p(2)) ~= 0
            error('There is something wrong.');
        end
    end
    disp('It is all right.');
end
