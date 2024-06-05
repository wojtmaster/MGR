function [v] = diff_eq_wiener_ARX(a, b, K, u_ucz, u_wer, y_ucz, y_wer, params, F_10, h_20, F_1_start, F_1_end, R, kk, tau, n)
    y_mod_ucz = zeros(1, kk/2);
    y_mod_wer = zeros(1, kk/2);
    y_mod_ucz(1:tau+2) = y_ucz(1:tau+2);
    y_mod_wer(1:tau+2) = y_wer(1:tau+2);

    v = zeros(1, kk/2);

    % Model ARX
    % Zbiór uczący
    for k = tau+3:kk/2
        v(k) = - a'*[y_ucz(k-1:-1:k-2)]' + (b(tau+1:tau+2)/K)'*[u_ucz(1,k-(tau+1):-1:k-(tau+2))]' + b(tau+1:tau+2)'*[u_ucz(2, k-1:-1:k-2)]';
        index = round((v(k)+F_10-F_1_start)*n / (F_1_end-F_1_start));
        y_mod_ucz(k) = find_value(params, v(k)+F_10, index, R) - h_20;
    end
    
    v = zeros(1, kk/2);
    % Zbiór weryfikujący
    for k = tau+3:kk/2
        v(k) = - a'*[y_wer(k-1:-1:k-2)]' + (b(tau+1:tau+2)/K)'*[u_wer(1,k-(tau+1):-1:k-(tau+2))]' + b(tau+1:tau+2)'*[u_wer(2, k-1:-1:k-2)]';
        index = round((v(k)+F_10-F_1_start)*n / (F_1_end-F_1_start));
        y_mod_wer(k) = find_value(params, v(k)+F_10, index, R) - h_20;
    end

    E_ucz = sum((y_ucz - y_mod_ucz).^2)/(kk/2);
    E_wer = sum((y_wer - y_mod_wer).^2)/(kk/2);
    fprintf('Model ARX \n');
    fprintf('E_ucz = %.3f \n', E_ucz);
    fprintf('E_wer = %.3f \n', E_wer);

    figure;
    hold on;
    stairs(0:kk/2-1, y_mod_ucz, 'b-', 'LineWidth', 1.2);
    stairs(0:kk/2-1, y_ucz, 'r-', 'LineWidth', 0.8);
    hold off;
    xlabel('k');
    ylabel('y_{ucz}(k)');
    plot_title = sprintf('Zbiór uczący - y_{ucz}(k) \n E_{ucz} = %.3f', E_ucz);
    title(plot_title);
    legend('y_{mod}', 'y');
    grid on;

    figure;
    hold on;
    stairs(0:kk/2-1, y_mod_wer, 'b-', 'LineWidth', 1.2);
    stairs(0:kk/2-1, y_wer, 'r-', 'LineWidth', 0.8);
    hold off;
    xlabel('k');
    ylabel('y_{wer}(k)');
    plot_title = sprintf('Zbiór weryfikujący - y_{wer}(k) \n E_{wer} = %.3f', E_wer);
    title(plot_title);
    legend('y_{mod}', 'y');
    grid on;
end