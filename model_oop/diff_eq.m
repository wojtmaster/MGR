function diff_eq(a, b, u, y, kk, tau, mode)

    y_mod = zeros(size(y));
    y_mod(1:tau+2) = y(1:tau+2);

    if strcmp(mode, 'ARX')
        for k = tau+3:kk
            y_mod(k) = - a*[y(k-1:-1:k-2)]' + b*[u(1, k-(tau+1):-1:k-(tau+2))]' + b*[u(2, k-1:-1:k-2)]';
        end
    else
        for k = tau+3:kk
            y_mod(k) = - a*[y_mod(k-1:-1:k-2)]' + b*[u(1, k-(tau+1):-1:k-(tau+2))]' + b*[u(2, k-1:-1:k-2)]';
        end
    end
    
    E = sum((y - y_mod).^2)/(kk);
    fprintf('Model %s \n', mode);
    fprintf('E = %.3f \n', E);

    figure;
    hold on;
    stairs(0:kk-1, y_mod, 'b-', 'LineWidth', 1.2);
    stairs(0:kk-1, y, 'r-', 'LineWidth', 0.8);
    hold off;
    xlabel('k');
    ylabel('y(k)');
    plot_title = sprintf('Model %s \n E = %.3f', mode, E);
    title(plot_title);
    legend('y_{mod}', 'y');
    grid on;
    % file_name = sprintf('../raport/pictures/arx_ucz.pdf');
    % exportgraphics (gcf, file_name);

end