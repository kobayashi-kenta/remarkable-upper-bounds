function A = create_coefficient_matrix(quadratic_form, variables)
    m = length(variables);
    A = sym(zeros(m));
    for p = 1:m
        df = diff(quadratic_form, variables(p));
        for q = 1:p
            A(p, q) = simplifyFraction(diff(df, variables(q)) / 2);
            if q<p
                A(q, p) = A(p, q);
            end
        end
    end
end
