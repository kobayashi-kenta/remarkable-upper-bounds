function Lemma14_1
    pair = [];

    a = sym('a'); b = sym('b');
    phi = b^2/(a^2 + b^2);

    f = diff(phi, a);
    g = - 2*a*b^2/(a^2 + b^2)^2;
    h = - 2/sym(3)/b ...
        + (2*(a - b)^4 + 2*a*b*(2*a - b)^2)/3/b/(a^2 + b^2)^2;
    pair = [pair [f; g] [f; h]];

    f = diff(phi, a, 2);
    g = 2*b^2*(3*a^2 - b^2)/(a^2 + b^2)^3;
    pair = [pair [f; g]];

    f = diff(phi, b);
    g = 2*a^2*b/(a^2 + b^2)^2;
    pair = [pair [f; g]];

    f = diff(phi, b, 2);
    g = 2*a^2*(a^2 - 3*b^2)/(a^2 + b^2)^3;
    h = - 7/sym(9)/b^2 + ( ...
            3*a^2*b^2*(a^2 + b^2) ...
            + (7*a^2 + 3*b^2)*(2*a^2 - b^2)^2 ...
            + b^2*(13*a^2 - 5*b^2)^2 ...
        )/36/b^2/(a^2 + b^2)^3;
    pair = [pair [f; g] [f; h]];

    for p = pair
        if simplifyFraction(p(1) - p(2)) ~= 0
            error('There is something wrong.');
        end
    end
    disp('It is all right.');
end
