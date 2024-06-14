function Y = find_value(params, u, index, R, s)
    a = zeros(size(R));
    b = zeros(size(R));
    c = zeros(size(R));

    Y = 0;
    weight = 0;

    if strcmp(s, 'linear')
        for i = 1:length(R)
            a(i) = params(2*i-1);
            b(i) = params(2*i);
        end

        for i = 1:length(R)
            if index >= length(R{i})
                Y = Y + R{i}(end)*(a(i) + b(i)*u);
                weight = weight + R{i}(end);
            elseif index <= 0
                Y = Y + R{i}(1)*(a(i) + b(i)*u);
                weight = weight + R{i}(1);
            else
                Y = Y + R{i}(index)*(a(i) + b(i)*u);
                weight = weight + R{i}(index);
            end
        end
        Y = Y / weight;
    else
        for i = 1:length(R)
            a(i) = params(3*i-2);
            b(i) = params(3*i-1);
            c(i) = params(3*i);
        end
        for i = 1:length(R)
            if index >= length(R{i})
                % Y = Y + R{i}(end)*(a(i)*tanh(b(i)*u/10 + c(i)));
                Y = Y + R{i}(end)*(a(i)*sinh(b(i)*u + c(i)));
                weight = weight + R{i}(end);
            elseif index <= 0
                % Y = Y + R{i}(1)*(a(i)*tanh(b(i)*u/10 + c(i)));
                Y = Y + R{i}(1)*(a(i)*sinh(b(i)*u + c(i)));
                weight = weight + R{i}(1);
            else
                % Y = Y + R{i}(index)*(a(i)*tanh(b(i)*u/10 + c(i)));
                Y = Y + R{i}(index)*(a(i)*sinh(b(i)*u + c(i)));
                weight = weight + R{i}(index);
            end
        end
        Y = Y / weight;
    end
end