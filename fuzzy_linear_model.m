function Y = fuzzy_linear_model(params, F_1, F_1_start, F_1_end, R, n)
    a = zeros(size(R));
    b = zeros(size(R));
    for i = 1:length(R)
        a(i) = params(2*i-1);
        b(i) = params(2*i);
    end

    Y = zeros(size(F_1));

    for i = 1:length(F_1)
        index = round((F_1(i)-F_1_start)*n / (F_1_end-F_1_start));
        if index == 0
            index = 1;
        end
        weight = 0;
        for j = 1:length(R)
            Y(i) = Y(i) + R{j}(index)*(a(j) + b(j)*F_1(i));
            weight = weight + R{j}(index);
        end
        Y(i) = Y(i) / weight;
    end
end