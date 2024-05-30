function [] = PID_controller_fuzzy(Kp, Ti, Td, a, b, u, y_zad, T, kk, tau, t, R_U, R_Y, optimal_params_U, optimal_params_Y)
    r0 = Kp * (1 + T/(2*Ti) + Td/T);
    r1 = Kp * (T/(2*Ti) - 2*Td/T - 1);
    r2 = Kp*Td/T;
    
    F_10 = 90;
    h_20 = 36;
    u_0 = 45;
    y_0 = 10;
    
    u(1,:) = zeros(1,kk);
    u_static = zeros(1,kk);
    y = zeros(1,kk);
    y_static = zeros(1,kk);
    e = zeros(kk,1);
    
    for k = 2:kk
        if k == 2
            y(k) = - a(1)*y(k-1) + b(1)*u(2,k-1);
        elseif k <= tau+1
            y(k) = - a(1)*y(k-1) - a(2)*y(k-2) + b(1)*u(2,k-1) + b(2)*u(2,k-2);
        elseif k == tau+2
            u_static(k-(tau+1)) = find_value(optimal_params_Y, y_static(k-1)+h_20, y_0, R_Y);
            u_static(k-(tau+1)) = u_static(k-(tau+1)) - F_10; 
            y(k) = - a(1)*y(k-1) - a(2)*y(k-2) + b(1)*u_static(k-(tau+1)) + b(1)*u(2,k-1) + b(2)*u(2,k-2);
        else
            u_static(k-(tau+1)) = find_value(optimal_params_Y, y_static(k-1)+h_20, y_0, R_Y);
            u_static(k-(tau+1)) = u_static(k-(tau+1)) - F_10; 
            y(k) = - a(1)*y(k-1) - a(2)*y(k-2) + b(1)*u_static(k-(tau+1)) + b(2)*u_static(k-(tau+2)) + b(1)*u(2,k-1) + b(2)*u(2,k-2);
        end
        
        % Obliczanie uchybu
        e(k) = y_zad(k) - y(k);
        
        % Obliczanie sterowania
        if k == 2
            u(1,k) = u(1,k-1) + r0*e(k) + r1*e(k-1);
        else
            u(1,k) = u(1,k-1) + r0*e(k) + r1*e(k-1) + r2*e(k-2);
        end

        y_static(k) = find_value(optimal_params_U, u(1,k)+F_10, u_0, R_U);
        y_static(k) = y_static(k) - h_20;
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