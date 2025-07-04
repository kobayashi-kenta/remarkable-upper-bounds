function z = L3ab(a, b)
    AA = (1-a)^2 + b^2;
    BB = a^2 + b^2;
    CC = 1;
    SS = b^2 / 4;
    z = (AA*BB + BB*CC + CC*AA) / 83 ...
        - (AA*BB*CC / (AA + BB + CC) + SS) / 24;
end
