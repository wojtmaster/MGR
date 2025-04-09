function y = myFun(x, params)
    A = params(1);
    B = params(2);
    C = params(3);

    y = A*x(1)^2 - B*x(2)^2 - C;
end