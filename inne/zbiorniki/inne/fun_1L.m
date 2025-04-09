function [x] = fun_1L(h_1, F_1, F_D, alfa_1, h_10)
    x = F_1 + F_D - alfa_1 * sqrt(h_10) - alfa_1 / (2*sqrt(h_10)) * (h_1-h_10);
end