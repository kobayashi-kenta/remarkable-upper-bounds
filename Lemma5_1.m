function Lemma5_1
    pair = [];

    a = sym('a'); b = sym('b'); h = sym('h');
    w1 = sym('w1'); w2 = sym('w2'); w3 = sym('w3');
    u0 = sym('u0');

    x = sym('x'); y = sym('y');
    t = sym('t'); t1 = sym('t1'); t2 = sym('t2');

    Pu = 1/b/h*( ...
            (b*(2*x - h) + 2*(1 - a)*y)*w1 ...
            - (b*(2*x - h) - 2*a*y)*w2 ...
            - (2*y - b*h)*w3 ...
        ) + 2*(3*(x^2 + y^2) - 2*(1 + a)*h*x - 2*b*h*y + a*h^2) ...
            /(1 - a + a^2 + b^2)/h^2*(w1 + w2 + w3 - 3*u0);
    Pux = diff(Pu, x); Puy = diff(Pu, y);

    W1 = int(subs(Pu, {x, y}, {h*(1 - (1 - a)*t), h*b*t}), t, 0, 1);
    W2 = int(subs(Pu, {x, y}, {h*a*t, h*b*t}), t, 0, 1);
    W3 = int(subs(Pu, {x, y}, {h*t, 0}), t, 0, 1);
    U0 = 2*int(int(subs(Pu, {x, y}, {h*(t1 + a*t2), h*b*t2}), ...
            t1, 0, 1 - t2), t2, 0, 1);
    pair = [pair [W1; w1] [W2; w2] [W3; w3] [U0; u0]];

    s = b*h^2/2;
    c1 = 2*((1 - a)^2 + b^2)/b; c2 = 2*(a^2 + b^2)/b; c3 = 2/b;

    v = subs(Pu^2, {x, y}, {h*(t1 + a*t2), h*b*t2});
    L2 = 2*s*int(int(v, t1, 0, 1 - t2), t2, 0, 1);
    l2 = F_alpha_0(s, c1, c2, c3, w1, w2, w3, u0);
    pair = [pair [L2; l2]];

    v = subs(Pux^2 + Puy^2, {x, y}, {h*(t1 + a*t2), h*b*t2});
    H1 = 2*s*int(int(v, t1, 0, 1 - t2), t2, 0, 1);
    h1 = F_alpha_1(s, c1, c2, c3, w1, w2, w3, u0);
    pair = [pair [H1; h1]];

    for p = pair
        if simplifyFraction(p(1) - p(2)) ~= 0
            error('There is something wrong.');
        end
    end
    disp('It is all right.');
end
