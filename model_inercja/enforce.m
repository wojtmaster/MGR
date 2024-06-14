function [u] = enforce(kk)
    u = zeros(2,kk);
    % u(1, :) = F_1
    % u(1, :) = 1;
    for i = 1:kk/250
        u(1, (i-1)*250+1 : i*250) = (rand()*5 - 2.5) * 10;
    end
    for i = 1:kk
        if(u(1,i) < -90)
            u(1,i) = -90;
        end
    end
    
    % u(2, :) = F_D
    u(2, 1 : 0.25*kk) = 0.2;
    u(2, 0.25*kk+1 : 0.5*kk) = -0.2;
    u(2, 0.5*kk+1 : 0.75*kk) = -0.2;
    u(2, 0.75*kk+1 : end) = 0.1;
    for i = 1:kk
        if(u(2,i) < -30)
            u(2,i) = -30;
        end
    end

    figure;
    stairs(0:kk-1, u(1,:), 'r-', 'LineWidth', 1.2);
    hold on;
    stairs(0:kk-1, u(2,:), 'b--', 'LineWidth', 1.2);
    hold off;
    xlabel('k');
    ylabel('u(k)');
    title('Wartości sygnałów sterujących u(k)');
    legend('F_1(k)', 'F_D(k)');
    grid on;
end