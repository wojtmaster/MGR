function [R, optimal_params] = static_characteristic_y_u(F_D0, alpha_2, u_start, u_end)
    u = linspace(u_start, u_end, (u_end-u_start)+1);
    d = F_D0;
    y = ((u+d) / alpha_2).^2;
    
    fig_u = figure;
    figure(fig_u);
    plot(u, y, 'b-');
    hold on;
    
    % Fuzzy sets
    R = cell(1,5);
    alpha = [-0.1, 0.1, 0.1, 0.1, 0.1];
    c = [45 70 90 110 135];
    
    for i = 1:length(R)
        for j = 1:length(u)
            if i == 1 || i == length(R)
                R{i}(j) = 1 / (1 + exp(-alpha(i)*(u(j)-c(i))));
            else
                R{i}(j) = 1 / (1 + exp(-alpha(i)*(u(j)-c(i-1)))) - 1 / (1 + exp(-alpha(i)*(u(j)-c(i+1))));
            end
        end
    end
    
    figure;
    for i = 1:length(R)
        plot(u, R{i});
        hold on;
    end
    grid on;
    xlabel('u');
    ylabel('y');
    title('Zbiory rozmyte y(u)');
    
    a = ones(size(R));
    b = ones(size(R));
     
    % MNK: Minimalizacja sumy kwadratów błędów
    fun = @(params) sum((y - fuzzy_linear_model(params, u, R)).^2);
    initial_params = [a(1), b(1)];
    for i = 2:length(R)
        initial_params = [initial_params, a(i), b(i)];
    end
    opt = optimset('MaxFunEvals', 10^6, 'MaxIter', 10^6);
    optimal_params = fminsearch(fun, initial_params, opt);
    
    % Prezentacja wyników
    Y = fuzzy_linear_model(optimal_params, u, R);
    figure(fig_u);
    plot(u, Y, 'ro');
    axis([45 135 10 70]);
    grid on;
    xlabel('u');
    ylabel('y');
    title('Charakterystyka statyczna y(u)');
    legend('y(u)', 'y(u) - fuzzy', 'Location', 'northwest');
end