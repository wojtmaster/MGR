function [x] = fun_2L(h_1, h_2, alfa_1, alfa_2, h_10, h_20)
    x = alfa_1 * sqrt(h_10) - alfa_2 * sqrt(h_20) + alfa_1 / (2*sqrt(h_10)) * (h_1-h_10) - alfa_2 / (2*sqrt(h_20)) * (h_2-h_20);
end