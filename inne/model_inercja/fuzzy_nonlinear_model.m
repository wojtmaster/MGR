function Y = fuzzy_nonlinear_model(params, u, u_start, u_end, R, n)
    a = zeros(size(R));
    b = zeros(size(R));
    c = zeros(size(R));
    for i = 1:length(R)
        a(i) = params(3*i-2);
        b(i) = params(3*i-1);
        c(i) = params(3*i);
    end

    Y = zeros(size(u));

    for i = 1:length(u)
        index = round((u(i)-u_start)*n / (u_end-u_start));
        if index <= 0
            index = 1;
        end
        weight = 0;
        for j = 1:length(R)
            % Y(i) = Y(i) + R{j}(index)*(a(j) + tanh(b(j)*u(i)/10));
            % Y(i) = Y(i) + R{j}(index)*(a(j)*tanh(b(j)*u(i)/10 + c(j)));

            % Y(i) = Y(i) + R{j}(index)*(a(j)+sinh(b(j)*u(i)));
            Y(i) = Y(i) + R{j}(index)*(a(j)*sinh(b(j)*u(i) + c(j)));
            weight = weight + R{j}(index);
        end
        Y(i) = Y(i) / weight;
    end
end