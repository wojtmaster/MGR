function Y = fuzzy_nonlinear_model(params, F_1, R)
    a = zeros(size(R));
    b = zeros(size(R));
    for i = 1:length(R)
        a(i) = params(2*i-1);
        b(i) = params(2*i);
    end

    Y = zeros(size(F_1));

    for i = 1:length(F_1)
        weight = 0;
        for j = 1:length(R)
            Y(i) = Y(i) + R{j}(i)*(a(j) + tanh(b(j)*F_1(i)));
            % Y(i) = Y(i) + R{j}(i)*sinh(a(j) + b(j)*u(i));
            weight = weight + R{j}(i);
        end
        Y(i) = Y(i) / weight;
    end
end