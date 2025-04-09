%% Static charactericstic
clear;
% Przykładowe dane wejściowe i wyjściowe
y = linspace(10, 70, 61);
d = 30;
u = 20*sqrt(y) - d;

figure(1);
plot(y, u, 'b-');
hold on;
xlabel('y');
ylabel('u');
title('Charakterystyka statyczna u(y)');

% Fuzzy sets
R = cell(1,5);
a = 0.1;
alpha = [-a, a, a, a, a];
c = [10 25 40 55 70];

for i = 1:length(R)
    for j = 1:length(y)
        if i == 1 || i == length(R)
            R{i}(j) = 1 / (1 + exp(-alpha(i)*(y(j)-c(i))));
        else
            R{i}(j) = 1 / (1 + exp(-alpha(i)*(y(j)-c(i-1)))) - 1 / (1 + exp(-alpha(i)*(y(j)-c(i+1))));
        end
    end
end

% figure;
% for i = 1:length(R)
%     plot(y, R{i});
%     hold on;
% end

% Optimization
a = ones(size(R));
b = ones(size(R));
 
% MNK: Minimalizacja sumy kwadratów błędów
fun = @(params) sum((u - fuzzy_linear_model(params, y, R)).^2);
initial_params = [a(1), b(1)];
for i = 2:length(R)
    initial_params = [initial_params, a(i), b(i)];
end
opt = optimset('MaxFunEvals', 10^6, 'MaxIter', 10^6);
optimal_params = fminsearch(fun, initial_params, opt);

% Prezentacja wyników
U = fuzzy_linear_model(optimal_params, y, R);
figure(1);
plot(y, U, 'ro');

% Znajdowanie wartości u na charakterystyce statycznej
y = 25;
y_0 = 10;
u = find(optimal_params, y, y_0, R);
disp(u);

% Funkcja modelu Takagi-Sugeno
function U = fuzzy_linear_model(params, y, R)
    for i = 1:length(R)
        a(i) = params(2*i-1);
        b(i) = params(2*i);
    end

    U = zeros(size(y));

    for i = 1:length(y)
        weight = 0;
        for j = 1:length(R)
            U(i) = U(i) + R{j}(i)*(a(j) + b(j)*y(i));
            weight = weight + R{j}(i);
        end
        U(i) = U(i) / weight;
    end
end

function U = find(params, y, y_0, R)
    for i = 1:length(R)
        a(i) = params(2*i-1);
        b(i) = params(2*i);
    end

    U = 0;
    weight = 0;
    for i = 1:length(R)
        U = U + R{i}(y-(y_0-1))*(a(i) + b(i)*y);
        weight = weight + R{i}(y-(y_0-1));
    end
    U = U / weight;
end

function U = fuzzy_nonlinear_model(params, y, R)
    for i = 1:length(R)
        a(i) = params(2*i-1);
        b(i) = params(2*i);
    end

    U = zeros(size(y));

    for i = 1:length(y)
        weight = 0;
        for j = 1:length(R)
            % Y(i) = Y(i) + R{j}(i)*sinh(a(j)+b(j)*u(i));
            U(i) = U(i) + R{j}(i)*a(j)*tanh(b(j)*y(i));
            weight = weight + R{j}(i);
        end
        U(i) = U(i) / weight;
    end
end