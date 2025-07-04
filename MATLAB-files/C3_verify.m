function C3_verify
    format long g;

    %%  Parameters  %%
    n = 20;
    ha = 1/sym(50); hb = 1/sym(50);
    ca = sym(6); cb = sym(8);

    disp('Step 1: Create index');

    [vertex_list, max_index] = create_vertex_list(n, 1);
    [edge_list, max_index] = create_edge_list(n, max_index + 1);

    %%  Variable "constraints" is an array of the index  %%
    %%  that constraints are imposed.                    %%
    constraints ...
    = [vertex_list(1, 1) vertex_list(n+1, 1) vertex_list(1, n+1)];

    %%  "triangle(i,k), k=1,2,3,4,5,6" denotes indexes      %%
    %%  of u1,u2,u3,w1,w2,w3 associate with i-th triangle.  %%
    [triangle, orientation] ...
    = set_triangle_beta(n, vertex_list, edge_list);

    disp('Step 2: Create coefficient matrix');

    a = sym('a'); b = sym('b'); h = sym('h');
    u1 = sym('u1'); u2 = sym('u2'); u3 = sym('u3');
    w1 = sym('w1'); w2 = sym('w2'); w3 = sym('w3');
    variables = [u1 u2 u3 w1 w2 w3];
    e = sym(ones(3, 3));

    s = b/2*h^2;
    c1 = 2*((1-a)^2+b^2)/b; c2 = 2*(a^2+b^2)/b; c3 = 2/b;

    z = F_beta_0(s, c1, c2, c3, u1, u2, u3, w1, w2, w3);
    A1 = create_coefficient_matrix(z, variables);
    A2 = A1.*[e -e; -e e];

    z = F_beta_2(s, c1, c2, c3, u1, u2, u3, w1, w2, w3);
    B1 = create_coefficient_matrix(z, variables);
    B2 = B1.*[e -e; -e e];

	lambda = L3ab(a, b) / sym(n^4) * sym(n^4-1) ...
	        / (1+ca*ha^2) / (1+cb*hb^2);
	M1ab = lambda * B1 - A1;
    M2ab = lambda * B2 - A2;

	lambda = L3ab(a, 0) / sym(n^4) * sym(n^4-1) / (1+ca*ha^2);
    M1a0 = lambda * B1 - A1;
    M2a0 = lambda * B2 - A2;

    disp('Step 3: Verify the inequality for 1/10 <= b <= 1');

    h = 1 / intval(n);

    yk = 1000;
    for k = 1:119
        b = 1000 / intval(yk);
        xk = ceil(yk/40);
        for l = 0:xk
            a = l / intval(2*xk);
	        M1 = eval(M1ab);
    	    M2 = eval(M2ab);
            verify();
        end
        yk = floor(yk*51/50);
    end

    disp('Step 4: Verify the inequality for 0 < b <= 1/10');

    b = 1 / intval(10);
    for l = 0:250
        a = l / intval(500);
        M1 = eval(M1a0);
   	    M2 = eval(M2a0);
  	    verify();
    end

    disp('Verificaion complited');

    %%  Nested function for the verification  %%
    function verify()
        disp(['a=', num2str(mid(a)), ', b=', num2str(mid(b))]);
        CM = intval(zeros(max_index));
        for i = 1:n^2
            c = triangle(i, :);
            if orientation(i) == 1
                CM(c, c) = CM(c, c) + M1;
            else
                CM(c, c) = CM(c, c) + M2;
            end
        end

        %%  Impose constraints  %%
        CM(constraints, :) = [];
        CM(:, constraints) = [];

        if isspd(CM) == 1
            disp('Verified');
        else
            error('Verification failed.');
        end
    end
end
