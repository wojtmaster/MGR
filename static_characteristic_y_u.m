function [R, optimal_params] = static_characteristic_y_u(F_1, h_2, n, s)
    % F_1 = linspace(F_1_start, F_1_end, n);
    % F_D = F_D0;
    % y = ((F_1+F_D) / alpha_2).^2;
    
    % Fuzzy sets
    R = cell(1,5);
    alpha = [-0.1, 0.1, 0.1, 0.1, 0.1];
    c = [45 70 90 110 135];
    
    for i = 1:length(R)
        for j = 1:length(F_1)
            if i == 1 || i == length(R)
                R{i}(j) = 1 / (1 + exp(-alpha(i)*(F_1(j)-c(i))));
            else
                R{i}(j) = 1 / (1 + exp(-alpha(i)*(F_1(j)-c(i-1)))) - 1 / (1 + exp(-alpha(i)*(F_1(j)-c(i+1))));
            end
        end
    end
    
    figure;
    for i = 1:length(R)
        plot(F_1, R{i});
        hold on;
    end
    grid on;
    xlabel('u');
    ylabel('y');
    xlim([45 135]);
    title('Zbiory rozmyte y(u)');
    
    a = ones(size(R));
    b = ones(size(R));
     
    fig_u = figure;
    figure(fig_u);
    plot(F_1, h_2, 'b-', 'LineWidth', 1.2);
    hold on;

    if strcmp(s, 'linear')
        % MNK: Minimalizacja sumy kwadratów błędów
        fun = @(params) sum((h_2 - fuzzy_linear_model(params, F_1, R)).^2);
        initial_params = [a(1), b(1)];
        for i = 2:length(R)
            initial_params = [initial_params, a(i), b(i)];
        end
        opt = optimset('MaxFunEvals', 10^6, 'MaxIter', 10^6);
        optimal_params = fminsearch(fun, initial_params, opt);
        
        % Prezentacja wyników
        H_2 = fuzzy_linear_model(optimal_params, F_1, R);
        figure(fig_u);
        plot(F_1, H_2, 'ro');
        xlim([45 135]);
        ylim([10 70]);
        grid on;
        xlabel('F_1');
        ylabel('h_2');
        plot_title = sprintf('Charakterystyka statyczna h_2(F_1) \n Następniki liniowe');
        title(plot_title);
        legend('h_2(F_1)', 'h_2(F_1) - fuzzy', 'Location', 'northwest');
    else
        % MNK: Minimalizacja sumy kwadratów błędów
        fun = @(params) sum((h_2 - fuzzy_nonlinear_model(params, F_1, R)).^2);
        initial_params = [a(1), b(1)];
        for i = 2:length(R)
            initial_params = [initial_params, a(i), b(i)];
        end
        opt = optimset('MaxFunEvals', 10^6, 'MaxIter', 10^6);
        optimal_params = fminsearch(fun, initial_params, opt);
        
        % Prezentacja wyników
        H_2 = fuzzy_nonlinear_model(optimal_params, F_1, R);
        figure(fig_u);
        plot(F_1, H_2, 'ro');
        xlim([45 135]);
        ylim([10 70]);
        grid on;
        xlabel('F_1');
        ylabel('h_2');
        plot_title = sprintf('Charakterystyka statyczna h_2(F_1) \n Następniki hiperboliczne');
        title(plot_title);
        legend('y(u)', 'y(u) - fuzzy', 'Location', 'northwest');
    end

    error = sum((h_2 - H_2).^2)/n;
    disp('MSE = ');
    disp(error);

end