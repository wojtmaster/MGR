function [x] = fun_1(h_1, F_1, F_D, alfa_1)
    x = F_1 + F_D - alfa_1 * sqrt(h_1);
end