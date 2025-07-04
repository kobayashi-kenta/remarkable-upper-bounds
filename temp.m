function Lemma14_6
    pair = [];

    g = subs(omega_bb, b, bt)/48;

    f = 1008*( ...
            2*(24*(1 - a + a^2 + b^2)^2 - 39*b^2 )/2016 ...
            + b*(- (1 - 2*a)*(28 + 13*b^2 + 35*a*(1 - a))/1008) ...
        );

end
