function Lemma14_8
    pair = [];

    a = sym('a'); b = sym('b');
    d1 = 1 - a; d2 = 1 - 2*a;

    f = 112*(a^2 + b^2)*(d1^2 + b^2)/b^2*(L1ab(a,b) - L1ab(a,0));
    g = 8*(a*d1 - b^2)^2 + b^2;
    pair = [pair [f; g]];

    f = 864*(a^2 + b^2)*(d1^2 + b^2)/b^2*(L2ab(a,b) - L2ab(a,0));
    g = 32*(a*d1 - b^2)^2 + 5*b^2;
    pair = [pair [f; g]];

    f = 7968*(d1 + a^2 + b^2)*(d1 + a^2)/b^2 ...
        *(L3ab(a,b) - L3ab(a,0));
    g = 52*a*d1 + 39*d2^2 + 3*a^2*d1^2*(125 + 16*d2^2) ...
         + 3*b^2*(d1 + a^2)*(21 + 32*b^2 + 24*d2^2);
    pair = [pair [f; g]];

    for p = pair
        if simplifyFraction(p(1) - p(2)) ~= 0
            error('There is something wrong.');
        end
    end
    disp('It is all right.');
end
