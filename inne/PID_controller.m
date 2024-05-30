function [] = PID_controller(Kp, Ti, Td, a, b, u, y_zad, T, kk, tau, t)
    r0 = Kp * (1 + T/(2*Ti) + Td/T);
    r1 = Kp * (T/(2*Ti) - 2*Td/T - 1);
    r2 = Kp*Td/T;
    
    u(1,:) = zeros(1,kk);
    y = zeros(1,kk);
    e = zeros(kk,1);
    
    for i = 2:kk
        if i == 2
            y(i) = - a(1)*y(i-1) + b(1)*u(2,i-1);
        elseif i <= tau+1
            y(i) = - a(1)*y(i-1) - a(2)*y(i-2) + b(1)*u(2,i-1) + b(2)*u(2,i-2);
        elseif i == tau+2
            y(i) = - a(1)*y(i-1) - a(2)*y(i-2) + b(1)*u(1,i-(tau+1)) + b(1)*u(2,i-1) + b(2)*u(2,i-2);
        else
            y(i) = - a(1)*y(i-1) - a(2)*y(i-2) + b(1)*u(1,i-(tau+1)) + b(2)*u(1,i-(tau+2)) + b(1)*u(2,i-1) + b(2)*u(2,i-2);
        end
        
        % Obliczanie uchybu
        e(i) = y_zad(i) - y(i);
        
        % Obliczanie sterowania
        if i == 2
            u(1,i) = u(1,i-1) + r0*e(i) + r1*e(i-1);
        else
            u(1,i) = u(1,i-1) + r0*e(i) + r1*e(i-1) + r2*e(i-2);
        end
    end
    
    figure;
    hold on;
    stairs(0:kk-1, y_zad);
    stairs(0:kk-1, y);
    hold off;
    title('y(k)');
    xlabel('k');
    ylabel('y(k)');
    legend('y_{zad}', 'y');
    
    figure;
    stairs(0:kk-1, u(1,:));
    title('u(k)');
    xlabel('k');
    ylabel('u(k)');
    legend('u');
end