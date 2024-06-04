function [] = DMC_SL(a, b, N, Nu, D, lambda, u, y_zad, y_max, u_max, delta_u_max, A, C, alpha_1, alpha_2, h_10, h_20, kk, tau, Tp)

    M = zeros(N, Nu);
    M_p = zeros(N, D-1);
    
    % Alokacja pamięci
    y = zeros(2, kk);
    error = zeros(1, kk);
    
    delta_up = zeros(D-1,1);
    delta_u = 0;
    u(1, :) = zeros(1, kk);
    
    % Sterowanie DMC
    for k = 2:kk
        if k == 2
            y(1,k) = - a(1,1)*y(1,k-1) + b(tau+1,1)*u(2,k-1);
            y(2,k) = - a(1,2)*y(2,k-1) + b(tau+1,2)*u(2,k-1);
            
        elseif k <= tau+1
            y(1,k) = - a(1,1)*y(1,k-1) + b(tau+1,1)*u(2,k-1);
            y(2,k) = - a(1,2)*y(2,k-1) - a(2,2)*y(2,k-2) + b(tau+1,2)*u(2,k-1) + b(tau+2,2)*u(2,k-2);
            
        elseif k == tau+2
            y(1,k) = - a(1,1)*y(1,k-1) + b(tau+1,1)*u(1,k-(tau+1)) + b(tau+1,1)*u(2,k-1);
            y(2,k) = - a(1,2)*y(2,k-1) - a(2,2)*y(2,k-2) + b(tau+1,2)*u(1,k-(tau+1)) + b(tau+1,2)*u(2,k-1) + b(tau+2,2)*u(2,k-2);
        else
            y(1,k) = - a(1,1)*y(1,k-1) + b(tau+1,1)*u(1,k-(tau+1)) + b(tau+1,1)*u(2,k-1);
            y(2,k) = - a(1,2)*y(2,k-1) - a(2,2)*y(2,k-2) + b(tau+1,2)*u(1,k-(tau+1)) + b(tau+2,2)*u(1,k-(tau+2)) + b(tau+1,2)*u(2,k-1) + b(tau+2,2)*u(2,k-2);
        end
    
        % Rzutowanie na ograniczenia
        if(y(1,k) > y_max)
            y(1,k) = y_max;
        elseif(y(1,k) < -y_max)
            y(1,k) = -y_max;
        end

        if(y(2,k) > y_max)
            y(2,k) = y_max;
        elseif(y(2,k) < -y_max)
            y(2,k) = -y_max;
        end

        % Tutaj trzeba przejść na zmienne F_1, F_D, h_1 i h_2
        % Obliczyć na nowo wszytskie współczynniki i póścić algorytm dalej
        V_10 = A * (y(1,k) + h_10);
        V_20 = C * (y(2,k) + h_20)^2;
        G_z = tf_function(A, C, alpha_1, alpha_2, V_10, V_20, tau*Tp, Tp);
        a(1,1) = [G_z.Denominator{1,1}(2)];
        b(tau+1,1) = [G_z.Numerator{1,1}(2)];
        % Współczynniki do równań różnicowych na h_2
        b(tau+1,2) = [G_z.Numerator{2,1}(2)];
        b(tau+2,2) = [G_z.Numerator{2,1}(3)];
        a(1,2) = [G_z.Denominator{2,1}(2)];
        a(2,2) = [G_z.Denominator{2,1}(3)];

        % % Odp. skokowa
        s = step(G_z(2,1), Tp*(D-1));

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

        % Przepisanie sterowań do wektora przeszłych sterowań
        delta_up = [delta_u; delta_up(1:end-1)];
    
        % Oblicznie uchybu    
        e = y_zad(k) - y(2,k);
    
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
        error(k) = (y_zad(k) - y(2,k))^2;
    end
    disp(sum(error)/kk);
    
    %% Prezentacja wyników
    figure;
    hold on;
    stairs(0:kk-1, y_zad, 'LineWidth', 1.2);
    stairs(0:kk-1, y(2,:), 'LineWidth', 1.2);
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