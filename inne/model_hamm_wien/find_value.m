function Y = find_value(params, u, index, R)
    a = zeros(size(R));
    b = zeros(size(R));
    for i = 1:length(R)
        a(i) = params(2*i-1);
        b(i) = params(2*i);
    end

    Y = 0;
    weight = 0;
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
end