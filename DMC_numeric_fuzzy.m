function [delta_u] = DMC_numeric_fuzzy(a, b, N, Nu, D, lambda, s, u, y_zad, y_max, u_max, delta_u_max, kk, tau, t, R_U, R_Y, optimal_params_U, optimal_params_Y)
    M = zeros(N, Nu);
    M_p = zeros(N, D-1);

    F_10 = 90;
    h_20 = 36;
    u_0 = 45;
    y_0 = 10;
    
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
    y_static = zeros(1, kk);
    error = zeros(1, kk);
    
    delta_up = zeros(D-1,1);
    delta_u = zeros(Nu,1);
    u(1, :) = zeros(1, kk);
    u_static = zeros(1, kk);
    
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
        %F_1
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
           
        % Przepisanie sterowań do wektora przeszłych sterowań
        for i = (D-1):-1:1
            if i == 1
                delta_up(i) = delta_u(1);
            else
                delta_up(i) = delta_up(i-1);
            end
        end
    
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
        
        % Obliczenie sterowania dla chwili (i+1) w chwili i-tej
        u(1,k) = u(1,k-1) + delta_u(1);

        % Odczyt wartości  
        y_static(k) = find_value(optimal_params_U, u(1,k)+F_10, u_0, R_U);
        y_static(k) = y_static(k) - h_20;
    
        % Obliczenie wskaźnika regulacji
        error(k) = error(k-1) + 1/2 * (y_zad(k) - y(k))^2;
    end
    disp(sum(error));
    
    % Prezentacja wyników
    figure;
    hold on;
    stairs(0:kk-1, y_zad);
    stairs(0:kk-1, y);
    hold off;
    legend('y_{zad}', 'y');
    title('Poziom wody w drugim zbiorniku - h_2');
    xlabel('k');
    ylabel('y(k)');
    
    figure;
    hold on;
    stairs(0:kk-1, u(1,:));
    hold off;
    legend('u');
    title('Sygnał sterujący F_1');
    xlabel('k');
    ylabel('u');
end