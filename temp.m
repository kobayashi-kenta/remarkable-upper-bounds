function L(a, b)
    AA = (1-a)^2 + b^2;
    BB = a^2 + b^2;
    CC = 1;
    SS = b^2 / 4;
    z = AA*BB*CC / 16 / SS - (AA + BB + CC) / 30 ...
        - SS / 5 * (1/AA + 1/BB + 1/CC);
end
