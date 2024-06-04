function [] = DMC_analitic(a, b, N, Nu, D, lambda, s, u, y_zad, y_max, u_max, delta_u_max, kk, tau)

    M = zeros(N, Nu);
    M_p = zeros(N, D-1);

    % Implementacja macierzy M_p
    for i = 1:N
        for j = 1:D-1
            if(i+j <= D)
                M_p(i,j) = s(i+j) - s(j);
            else 
                M_p(i,j) = s(D) - s(j);
            end
        end
    end
    
    % Implementacja macierzy M
    for i = 1:N
        for j = 1:Nu
            if(i >= j)
                M(i,j) = s(i-j+1);
            end
        end
    end

    % Macierz K
    K = ((M' * M + lambda * eye(Nu))^(-1)) * M';
    ke = sum(K(1, :));
    ku = K(1, :) * M_p;
    
    % Alokacja pamięci
    y = zeros(1, kk);
    error = zeros(1, kk);
    
    delta_up = zeros(D-1,1);
    delta_u = 0;
    u(1, :) = zeros(1, kk);
    
    % Sterowanie DMC
    for k = 2:kk
        if k == 2
            y(k) = - a(1)*y(k-1) + b(tau+1)*u(2,k-1);
        elseif k <= tau+1
            y(k) = - a(1)*y(k-1) - a(2)*y(k-2) + b(tau+1)*u(2,k-1) + b(tau+2)*u(2,k-2);
        elseif k == tau+2
            y(k) = - a(1)*y(k-1) - a(2)*y(k-2) + b(tau+1)*u(1,k-(tau+1)) + b(tau+1)*u(2,k-1) + b(tau+2)*u(2,k-2);
        else
            y(k) = - a(1)*y(k-1) - a(2)*y(k-2) + b(tau+1)*u(1,k-(tau+1)) + b(tau+2)*u(1,k-(tau+2)) + b(tau+1)*u(2,k-1) + b(tau+2)*u(2,k-2);
        end
    
        % Rzutowanie na ograniczenia
        if(y(k) > y_max)
            y(k) = y_max;
        elseif(y(k) < -y_max)
            y(k) = -y_max;
        end
    
        % Przepisanie sterowań do wektora przeszłych sterowań
        delta_up = [delta_u; delta_up(1:end-1)];
    
        % Oblicznie uchybu    
        e = y_zad(k) - y(k);
    
        % Obliczenie przyrostu sterowania dla chwili (i+1) w chwili i-tej
        delta_u = ke * e - ku * delta_up;
        
        % Ograniczenie wartości przyrostu sterowania
        if(delta_u > delta_u_max)
            delta_u = delta_u_max;
        elseif(delta_u < -delta_u_max)
            delta_u = -delta_u_max;
        end
        
        % Obliczenie sterowania
        u(1,k) = u(1,k-1) + delta_u; 
        
        % Ograniczenie sterowania
        if(u(1,k) > u_max)
            u(1,k) = u_max;
        elseif(u(1,k) < -u_max)
            u(1,k) = -u_max;
        end
    
        % Obliczenie wskaźnika regulacji
        error(k) = (y_zad(k) - y(k))^2;
    end
    disp(sum(error)/kk);
    
    %% Prezentacja wyników
    figure;
    hold on;
    stairs(0:kk-1, y_zad, 'LineWidth', 1.2);
    stairs(0:kk-1, y, 'LineWidth', 1.2);
    hold off;
    legend('y_{zad}', 'y');
    title('Poziom wody w drugim zbiorniku - h_2');
    xlabel('k');
    ylabel('y(k)');
    grid on;
    
    figure;
    hold on;
    stairs(0:kk-1, u(1,:), 'LineWidth', 1.2);
    hold off;
    legend('u');
    title('Sygnał sterujący F_1');
    xlabel('k');
    ylabel('u');
    grid on;
end