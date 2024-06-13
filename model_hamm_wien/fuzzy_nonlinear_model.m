function Y = fuzzy_nonlinear_model(params, u, u_start, u_end, R, n)
    a = zeros(size(R));
    b = zeros(size(R));
    for i = 1:length(R)
        a(i) = params(2*i-1);
        b(i) = params(2*i);
    end

    Y = zeros(size(u));

    for i = 1:length(u)
        index = round(u(i)*n / (u_end-u_start));
        if index == 0
            index = 1;
        end
        weight = 0;
        for j = 1:length(R)
            % Y(i) = Y(i) + R{j}(index)*(a(j) + tanh(b(j)*F_1(i)));
            % Y(i) = Y(i) + R{j}(index)*(a(j)+sinh(b(j)*F_1(i)));
            Y(i) = Y(i) + R{j}(index)*(a(j)+b(j)*u(i)^2);
            weight = weight + R{j}(index);
        end
        Y(i) = Y(i) / weight;
    end
end