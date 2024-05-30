function Y = fuzzy_linear_model(params, u, R)
    a = zeros(size(R));
    b = zeros(size(R));
    for i = 1:length(R)
        a(i) = params(2*i-1);
        b(i) = params(2*i);
    end

    Y = zeros(size(u));

    for i = 1:length(u)
        weight = 0;
        for j = 1:length(R)
            Y(i) = Y(i) + R{j}(i)*(a(j) + b(j)*u(i));
            weight = weight + R{j}(i);
        end
        Y(i) = Y(i) / weight;
    end
end