function [delta_u] = DMC_numeric(a, b, N, Nu, D, lambda, s, u, y_zad, y_max, u_max, delta_u_max, kk, tau)
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
    
    % Alokacja pamięci
    y = zeros(1, kk);
    error = zeros(1, kk);
    
    delta_up = zeros(D-1,1);
    delta_u = zeros(Nu,1);
    u(1, :) = zeros(1, kk);
    
    Y_max = ones(N,1)*y_max;
    Y_min = ones(N,1)*(-y_max);
    
    U_max = ones(Nu,1)*u_max;
    U_min = ones(Nu,1)*(-u_max);
    delta_U_max = ones(Nu,1)*delta_u_max;
    delta_U_min = ones(Nu,1)*(-delta_u_max);
    
    J = tril(ones(Nu));
    A = [-J; J; -M; M];
    H = 2*(M'*M + lambda*eye(Nu));
    
    % Sterowanie DMC
    for k = 2:kk
        % % F_1
        if k == 2
            y(k) = - a(1)*y(k-1) + b(1)*u(2,k-1);
        elseif k <= tau+1
            y(k) = - a(1)*y(k-1) - a(2)*y(k-2) + b(1)*u(2,k-1) + b(2)*u(2,k-2);
        elseif k == tau+2
            y(k) = - a(1)*y(k-1) - a(2)*y(k-2) + b(1)*u(1,k-(tau+1)) + b(1)*u(2,k-1) + b(2)*u(2,k-2);
        else
            y(k) = - a(1)*y(k-1) - a(2)*y(k-2) + b(1)*u(1,k-(tau+1)) + b(2)*u(1,k-(tau+2)) + b(1)*u(2,k-1) + b(2)*u(2,k-2);
        end
           
        % Przepisanie sterowań do wektora przeszłych sterowań
        delta_up = [delta_u(1); delta_up(1:end-1)];
    
        U_k_1 = ones(Nu,1)*u(1,k-1);
        Y = ones(N,1)*y(k);
        Y_0 = Y + M_p*delta_up;
        Y_zad = ones(N,1)*y_zad(k);
        B = [-U_min + U_k_1;
            U_max - U_k_1;
            -Y_min + Y_0;
            Y_max - Y_0];
        f = -2*M'*(Y_zad - Y_0);
        
        % Obliczenie przyrostu sterowania dla chwili (i+1) w chwili i-tej
        delta_u = quadprog(H,f,A,B,[],[],delta_U_min,delta_U_max);
        
        % Obliczenie sterowania
        u(1,k) = u(1,k-1) + delta_u(1); 
    
        % Obliczenie wskaźnika regulacji
        error(k) = (y_zad(k) - y(k))^2;
    end
    disp(sum(error) / kk);
    
    % Prezentacja wyników
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