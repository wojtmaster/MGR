function [R, optimal_params] = static_characteristic_y_v(v, v_ucz, v_wer, y_ucz, y_wer, v_start, v_end, n, s)
    
    % Fuzzy sets
    if strcmp(s, 'linear')
        R = cell(1,5);
        sigma = [12 10 10 10 12];
        % mean = [-20 -10 0 10 20];
        mean = [-30 -15 0 15 30];
        for i = 1:length(R)
            R{i} = gaussmf(v, [sigma(i) mean(i)]);
        end
    else
        R = cell(1,3);
        sigma = 10;
        mean = [-20 0 20];
        for i = 1:length(R)
            R{i} = gaussmf(v, [sigma mean(i)]);
        end
    end
    
    figure;
    for i = 1:length(R)
        plot(v, R{i});
        hold on;
    end
    grid on;
    xlabel('v');
    ylabel('y');
    xlim([v_start v_end]);
    title('Zbiory rozmyte y(v)');
    
    a = ones(size(R));
    b = ones(size(R));
     
    fig_u_ucz = figure;
    figure(fig_u_ucz);
    plot(v_ucz, y_ucz, 'b-', 'LineWidth', 1.2);
    hold on;

    fig_u_wer = figure;
    figure(fig_u_wer);
    plot(v_wer, y_wer, 'b-', 'LineWidth', 1.2);
    hold on;

    if strcmp(s, 'linear')
        % MNK: Minimalizacja sumy kwadratów błędów
        fun = @(params) sum((y_ucz - fuzzy_linear_model(params, v_ucz, v_start, v_end, R, n)).^2);
        initial_params = [a(1), b(1)];
        for i = 2:length(R)
            initial_params = [initial_params, a(i), b(i)];
        end
        opt = optimset('MaxFunEvals', 10^6, 'MaxIter', 10^6);
        optimal_params = fminsearch(fun, initial_params, opt);
        
        % Prezentacja wyników
        % Zbiór uczący
        y_mod_ucz = fuzzy_linear_model(optimal_params, v_ucz, v_start, v_end, R, n);
        figure(fig_u_ucz);
        plot(v_ucz, y_mod_ucz, 'ro');
        xlim([v_start v_end]);
        ylim([-40 40]);
        grid on;
        xlabel('u');
        ylabel('y');
        plot_title = sprintf('Charakterystyka statyczna y_{ucz}(v_{ucz}) \n Następniki liniowe');
        title(plot_title);
        legend('y(v)', 'y(v) - fuzzy', 'Location', 'northwest');
    
        % Zbiór weryfikujący
        y_mod_wer = fuzzy_linear_model(optimal_params, v_wer, v_start, v_end, R, n);
        figure(fig_u_wer);
        plot(v_wer, y_mod_wer, 'ro');
        xlim([v_start v_end]);
        ylim([-40 40]);
        grid on;
        xlabel('v');
        ylabel('y');
        plot_title = sprintf('Charakterystyka statyczna y_{wer}(v_{wer}) \n Następniki liniowe');
        title(plot_title);
        legend('y(v)', 'y(v) - fuzzy', 'Location', 'northwest');
    else
        % MNK: Minimalizacja sumy kwadratów błędów
        fun = @(params) sum((y_ucz - fuzzy_nonlinear_model(params, v_ucz, v_start, v_end, R, n)).^2);
        initial_params = [a(1), b(1)];
        for i = 2:length(R)
            initial_params = [initial_params, a(i), b(i)];
        end
        opt = optimset('MaxFunEvals', 10^6, 'MaxIter', 10^6);
        optimal_params = fminsearch(fun, initial_params, opt);
        
        % Prezentacja wyników
        y_mod_ucz = fuzzy_nonlinear_model(optimal_params, v_ucz, v_start, v_end, R, n);
        figure(fig_u_ucz);
        plot(v_ucz, y_mod_ucz, 'ro');
        xlim([-v_start v_end]);
        ylim([-40 40]);
        grid on;
        xlabel('v');
        ylabel('y');
        plot_title = sprintf('Charakterystyka statyczna y_{ucz}(v_{ucz}) \n Następniki hiperboliczne');
        title(plot_title);
        legend('y(v)', 'y(v) - fuzzy', 'Location', 'northwest');

        y_mod_wer = fuzzy_nonlinear_model(optimal_params, v_wer, v_start, v_end, R, n);
        figure(fig_u_wer);
        plot(v_wer, y_mod_wer, 'ro');
        xlim([-v_start v_end]);
        ylim([-40 40]);
        grid on;
        xlabel('v');
        ylabel('y');
        plot_title = sprintf('Charakterystyka statyczna y_{wer}(v_{wer}) \n Następniki hiperboliczne');
        title(plot_title);
        legend('y(v)', 'y(v) - fuzzy', 'Location', 'northwest');
    end

    E_ucz = sum((y_ucz - y_mod_ucz).^2)/(n/2);
    E_wer = sum((y_wer - y_mod_wer).^2)/(n/2);
    
    fprintf('Wartość E_ucz = %.3f \n', E_ucz);
    fprintf('Wartość E_wer = %.3f \n', E_wer);

end