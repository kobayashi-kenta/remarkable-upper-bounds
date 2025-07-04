function Lemma14_2
    pair = [];

    a = sym('a'); b = sym('b');
    phi = b^2/(a^2 + b^2);
    psi = b^2*(1 + 2*a)/(1 + 4*b^2)*phi;

    f = diff(psi, a);
    g = 2*b^2/(1 + 4*b^2)*phi ...
        + b^2*(1 + 2*a)/(1 + 4*b^2)*diff(phi, a);
    pair = [pair [f; g]];

    f = diff(psi, a, 2);
    g = 4*b^2/(1 + 4*b^2)*diff(phi, a) ...
        + b^2*(1 + 2*a)/(1 + 4*b^2)*diff(phi, a, 2);
    pair = [pair [f; g]];

    f = 4*b^2/(1 + 4*b^2)*(-2)/3/b ...
        + b^2*(1 + 2*a)/(1 + 4*b^2)*(-2)/b^2;
    g = -2*(3 + 6*a + 4*b)/3/(1 + 4*b^2);
    pair = [pair [f; g]];

    f = diff(psi, b);
    g = 2*b*(1 + 2*a)/(1 + 4*b^2)^2*phi ...
        + b^2*(1 + 2*a)/(1 + 4*b^2)*diff(phi, b);
    pair = [pair [f; g]];

    f = 2*b*(1 + 2*a)/(1 + 4*b^2)^2 ...
        + b^2*(1 + 2*a)/(1 + 4*b^2)*1/2/b;
    g = (1 + 2*a)*(2 + 16*b - 9*b^2)/10/(1 + 4*b^2) ...
        - (1 + 2*a)*(1 - b)/50/(1 + 4*b^2)^2 ...
            *(11*b^3 + 10*(1 - 3*b)^2 + b*(5 - 13*b)^2);
    pair = [pair [f; g]];

    f = diff(psi, b, 2);
    g = 2*(1 + 2*a)*(1 - 12*b^2)/(1 + 4*b^2)^3*phi ...
        + 4*b*(1 + 2*a)/(1 + 4*b^2)^2*diff(phi, b) ...
        + b^2*(1 + 2*a)/(1 + 4*b^2)*diff(phi, b, 2);
    pair = [pair [f; g]];

    f = 2*(1 + 2*a)*(1 - 12*b^2 - (5/sym(4) - 3*b^2)^2) ...
            /(1 + 4*b^2)^3*phi ...
        - b^2*(1 + 2*a)/(1 + 4*b^2)*7/sym(9)/b^2;
    g = -9*(1 + 2*a)/8/(1 + 4*b^2)*phi ...
        - 7*(1 + 2*a)/9/(1 + 4*b^2);
    pair = [pair [f; g]];

    for p = pair
        if simplifyFraction(p(1) - p(2)) ~= 0
            error('There is something wrong.');
        end
    end
    disp('It is all right.');
end
