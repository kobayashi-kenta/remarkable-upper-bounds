function z = L1ab(a, b)
    AA = (1-a)^2 + b^2;
    BB = a^2 + b^2;
    CC = 1;
    SS = b^2 / 4;
    z = (AA + BB + CC) / 28 - SS^2 / (AA * BB * CC);
end
