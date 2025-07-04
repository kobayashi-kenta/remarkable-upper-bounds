function Lemma5_2
    pair = [];

    a = sym('a'); b = sym('b'); h = sym('h');
    u1 = sym('u1'); u2 = sym('u2'); u3 = sym('u3');
    w1 = sym('w1'); w2 = sym('w2'); w3 = sym('w3');

    x = sym('x'); y = sym('y');
    t = sym('t'); t1 = sym('t1'); t2 = sym('t2');

    Pu = 1/b/h*( ...
            - (b*(x - h) + (1 - a)*y)*u1 ...
            + (b*x - a*y)*u2 + y*u3 ...
        ) + (b*(x - h) + (1 - a)*y)*(b*x + (1 - a)*y) ...
            /((1 - a)^2 + b^2)/b^2/h^2 ...
            *(b*w1 + ((1 - a)^2 + b^2)*u1 ...
                + (a*(1 - a) - b^2)*u2 - (1 - a)*u3) ...
        + (b*(x - h) - a*y)*(b*x - a*y) ...
            /(a^2 + b^2)/b^2/h^2 ...
            * (b*w2 + (a*(1-a) - b^2)*u1 ...
                + (a^2 + b^2)*u2 - a*u3) ...
        + (y - b*h)*y/b^2/h^2 ...
            * (b*w3 - (1-a)*u1 - a*u2 + u3);
    Pux = diff(Pu, x); Puy = diff(Pu, y);
    Puxx = diff(Pux, x); Puxy = diff(Pux, y); Puyy = diff(Puy, y);

    U1 = subs(Pu, {x, y}, {0, 0});
    U2 = subs(Pu, {x, y}, {h, 0});
    U3 = subs(Pu, {x, y}, {h*a, h*b});
    pair = [pair [U1; u1] [U2; u2] [U3; u3]];

    nx = b*h; ny = (1-a)*h;
    W1 = int(subs(Pux*nx + Puy*ny, ...
            {x, y}, {h*(1-(1-a)*t), h*b*t}), t, 0, 1);
    nx = -b*h; ny = a*h;
    W2 = int(subs(Pux*nx + Puy*ny, {x, y}, {h*a*t, h*b*t}), t, 0, 1);
    nx = 0; ny = -h;
    W3 = int(subs(Pux*nx + Puy*ny, {x, y}, {h*t, 0}), t, 0, 1);
    pair = [pair [W1; w1] [W2; w2] [W3; w3]];

    s = b*h^2/2;
    c1 = 2*((1 - a)^2 + b^2)/b; c2 = 2*(a^2 + b^2)/b; c3 = 2/b;

    v = subs(Pu^2, {x, y}, {h*(t1 + a*t2), h*b*t2});
    L2 = 2*s*int(int(v, t1, 0, 1 - t2), t2, 0, 1);
    l2 = F_beta_0(s, c1, c2, c3, u1, u2, u3, w1, w2, w3);
    pair = [pair [L2; l2]];

    v = subs(Pux^2 + Puy^2, {x, y}, {h*(t1 + a*t2), h*b*t2});
    H1 = 2*s*int(int(v, t1, 0, 1 - t2), t2, 0, 1);
    h1 = F_beta_1(s, c1, c2, c3, u1, u2, u3, w1, w2, w3);
    pair = [pair [H1; h1]];

    v = subs(Puxx^2 + 2*Puxy^2 + Puyy^2, ...
            {x, y}, {h*(t1 + a*t2), h*b*t2});
    H2 = 2*s*int(int(v, t1, 0, 1 - t2), t2, 0, 1);
    h2 = F_beta_2(s, c1, c2, c3, u1, u2, u3, w1, w2, w3);
    pair = [pair [H2; h2]];

    for p = pair
        if simplifyFraction(p(1) - p(2)) ~= 0
            error('There is something wrong.');
        end
    end
    disp('It is all right.');
end
