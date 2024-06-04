function [R, optimal_params] = static_characteristic_y_u(F_1, F_1_ucz, F_1_wer, h_2_ucz, h_2_wer, F_1_start, F_1_end, n, s)
    
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
     
    fig_u_ucz = figure;
    figure(fig_u_ucz);
    plot(F_1_ucz, h_2_ucz, 'b-', 'LineWidth', 1.2);
    hold on;

    fig_u_wer = figure;
    figure(fig_u_wer);
    plot(F_1_wer, h_2_wer, 'b-', 'LineWidth', 1.2);
    hold on;

    if strcmp(s, 'linear')
        % MNK: Minimalizacja sumy kwadratów błędów
        fun = @(params) sum((h_2_ucz - fuzzy_linear_model(params, F_1_ucz, F_1_start, F_1_end, R, n)).^2);
        initial_params = [a(1), b(1)];
        for i = 2:length(R)
            initial_params = [initial_params, a(i), b(i)];
        end
        opt = optimset('MaxFunEvals', 10^6, 'MaxIter', 10^6);
        optimal_params = fminsearch(fun, initial_params, opt);
        
        % Prezentacja wyników
        % Zbiór uczący
        h_2_mod_ucz = fuzzy_linear_model(optimal_params, F_1_ucz, F_1_start, F_1_end, R, n);
        figure(fig_u_ucz);
        plot(F_1_ucz, h_2_mod_ucz, 'ro');
        xlim([45 135]);
        ylim([10 70]);
        grid on;
        xlabel('F_1');
        ylabel('h_2');
        plot_title = sprintf('Charakterystyka statyczna h_{2ucz}(F_{1ucz}) \n Następniki liniowe');
        title(plot_title);
        legend('h_2(F_1)', 'h_2(F_1) - fuzzy', 'Location', 'northwest');
    
        % Zbiór weryfikujący
        h_2_mod_wer = fuzzy_linear_model(optimal_params, F_1_wer, F_1_start, F_1_end, R, n);
        figure(fig_u_wer);
        plot(F_1_wer, h_2_mod_wer, 'ro');
        xlim([45 135]);
        ylim([10 70]);
        grid on;
        xlabel('F_1');
        ylabel('h_2');
        plot_title = sprintf('Charakterystyka statyczna h_{2wer}(F_{1wer}) \n Następniki liniowe');
        title(plot_title);
        legend('h_2(F_1)', 'h_2(F_1) - fuzzy', 'Location', 'northwest');
    else
        % MNK: Minimalizacja sumy kwadratów błędów
        fun = @(params) sum((h_2_ucz - fuzzy_nonlinear_model(params, F_1_ucz, F_1_start, F_1_end, R, n)).^2);
        initial_params = [a(1), b(1)];
        for i = 2:length(R)
            initial_params = [initial_params, a(i), b(i)];
        end
        opt = optimset('MaxFunEvals', 10^6, 'MaxIter', 10^6);
        optimal_params = fminsearch(fun, initial_params, opt);
        
        % Prezentacja wyników
        h_2_mod_ucz = fuzzy_nonlinear_model(optimal_params, F_1_ucz, F_1_start, F_1_end, R, n);
        figure(fig_u_ucz);
        plot(F_1_ucz, h_2_mod_ucz, 'ro');
        xlim([45 135]);
        ylim([10 70]);
        grid on;
        xlabel('F_1');
        ylabel('h_2');
        plot_title = sprintf('Charakterystyka statyczna h_{2ucz}(F_{1ucz}) \n Następniki hiperboliczne');
        title(plot_title);
        legend('y(u)', 'y(u) - fuzzy', 'Location', 'northwest');

        h_2_mod_wer = fuzzy_nonlinear_model(optimal_params, F_1_wer, F_1_start, F_1_end, R, n);
        figure(fig_u_wer);
        plot(F_1_wer, h_2_mod_wer, 'ro');
        xlim([45 135]);
        ylim([10 70]);
        grid on;
        xlabel('F_1');
        ylabel('h_2');
        plot_title = sprintf('Charakterystyka statyczna h_{2wer}(F_{1wer}) \n Następniki hiperboliczne');
        title(plot_title);
        legend('y(u)', 'y(u) - fuzzy', 'Location', 'northwest');
    end

    E_ucz = sum((h_2_ucz - h_2_mod_ucz).^2)/(n/2);
    E_wer = sum((h_2_wer - h_2_mod_wer).^2)/(n/2);
    
    fprintf('Wartość E_ucz = %.3f \n', E_ucz);
    fprintf('Wartość E_wer = %.3f \n', E_wer);

end