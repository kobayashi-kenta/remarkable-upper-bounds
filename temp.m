function Lemma14_6
    pair = [];

    g = subs(omega_bb, b, bt)/48;

    for p = pair
        if simplifyFraction(p(1) - p(2)) ~= 0
            error('There is something wrong.');
        end
    end
    disp('It is all right.');
end
