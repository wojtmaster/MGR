function [R, optimal_params] = static_characteristic_u_y(F_D0, alpha_2, y_start, y_end)
    % Przykładowe dane wejściowe i wyjściowe
    y = linspace(y_start, y_end, (y_end-y_start)+1);
    d = F_D0;
    u = alpha_2*sqrt(y) - d;
    
    fig_y = figure;
    figure(fig_y);
    plot(y, u, 'b-');
    hold on;
    
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
    
    figure;
    for i = 1:length(R)
        plot(y, R{i});
        hold on;
    end
    grid on;
    xlabel('y');
    ylabel('u');
    title('Zbiory rozmyte u(y)');

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
    figure(fig_y);
    plot(y, U, 'ro');
    grid on;
    xlabel('y');
    ylabel('u');
    title('Charakterystyka statyczna u(y)');
    legend('u(y)', 'u(y) - fuzzy', 'Location', 'northwest');
    
    % Znajdowanie wartości u na charakterystyce statycznej
    % y = 25;
    % y_0 = 10;
    % u = find(optimal_params, y, y_0, R);
    % disp(u);
end