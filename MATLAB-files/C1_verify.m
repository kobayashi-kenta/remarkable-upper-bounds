function C1_verify
    format long g;

    %%  Parameters  %%
    n = 20;
    ha = 1 / sym(50); hb = 1 / sym(50);
    ca = sym(5); cb = sym(3);

    disp('Step 1: Create index');

    [edge_list, max_index] = create_edge_list(n, 1);
    [face_list, max_index] = create_face_list(n, max_index + 1);

    %%  "constraints" is an array of the index  %%
    %%  that constraints are imposed.           %%
    constraints = face_list(face_list > 0);

    %%  "triangle(i,k), k=1,2,3,4" denotes indexes    %%
    %%  of w1,w2,w3,u0 associate with i-th triangle.  %%
    triangle = set_triangle_alpha(n, edge_list, face_list);

    disp('Step 2: Create coefficient matrix');

    a = sym('a'); b = sym('b'); h = sym('h');
    w1 = sym('w1'); w2 = sym('w2'); w3 = sym('w3');
    u0 = sym('u0');
    variables = [w1 w2 w3 u0];

    s = b/2*h^2;
    c1 = 2*((1-a)^2+b^2)/b; c2 = 2*(a^2+b^2)/b; c3 = 2/b;

    z = F_alpha_0(s, c1, c2, c3, w1, w2, w3, u0);
    A = create_coefficient_matrix(z, variables);

    z = F_alpha_1(s, c1, c2, c3, w1, w2, w3, u0);
    B = create_coefficient_matrix(z, variables);

	lambda = L1ab(a, b) / sym(n^2) * sym(n^2-1) ...
	        / (1+ca*ha^2) / (1+cb*hb^2);
    Mab = lambda * B - A;

	lambda = L1ab(a, 0) / sym(n^2) * sym(n^2-1) / (1+ca*ha^2);
    Ma0 = lambda * B - A;

    disp('Step 3: Verify the inequality for 1/10 <= b <= 1');

    h = 1 / intval(n);

    yk = 1000;
    for k = 1:119
        b = 1000 / intval(yk);
        xk = ceil(yk/40);
        for l = 0:xk
            a = l / intval(2*xk);
            M = eval(Mab);
            verify();
        end
        yk = floor(yk*51/50);
    end

    disp('Step 4: Verify the inequality for 0 < b <= 1/10');

    b = 1 / intval(10);
    for l = 0:250
        a = l / intval(500);
        M = eval(Ma0);
  	    verify();
    end

    disp('Verificaion complited');

    %%  Nested function for the verification  %%
    function verify()
        disp(['a=', num2str(mid(a)), ', b=', num2str(mid(b))]);
        CM = intval(zeros(max_index));
        for i = 1:n^2
            c = triangle(i, :);
            CM(c, c) = CM(c, c) + M;
        end

        %%  Impose constraints  %%
        s = constraints(1);
        t = constraints(2:end);
        len = length(t);
        p = CM(:, s); q = CM(s, :); r = CM(s, s);
        CM(:, t) = CM(:, t) - repmat(p, 1, len);
        CM(t, :) = CM(t, :) - repmat(q, len, 1);
        CM(t, t) = CM(t, t) + repmat(r, len, len);
        CM(:, s) = []; CM(s, :) = [];

        if isspd(CM) == 1
            disp('Verified');
        else
            error('Verification failed.');
        end
    end
end
