function [y, y_L] = rk4(obj, u, kk)

    %% Alokacja pamięci
    y = zeros(1, kk);
    y_L = zeros(1, kk);
    
    %% Warunki początkowe
    V_1 = obj.V_10;
    V_2 = obj.V_20;
    V_1L = obj.V_10;
    V_2L = obj.V_20;
    
    %% Funkcje
    fun_1 = @(F_1, F_D, V_1) F_1 + F_D - obj.alpha_1 * (V_1/obj.A)^(1/2);
    fun_1L = @(F_1, F_D, V_1) F_1 + F_D - obj.alpha_1 * (obj.V_10/obj.A)^(1/2) - obj.alpha_1 / (2*obj.V_10^(1/2)*obj.A^(1/2)) * (V_1-obj.V_10);
    fun_2 = @(V_1, V_2) obj.alpha_1 * (V_1/obj.A)^(1/2) - obj.alpha_2 * (V_2/obj.C)^(1/4); 
    fun_2L = @(V_1, V_2) obj.alpha_1 * (obj.V_10/obj.A)^(1/2) - obj.alpha_2 * (obj.V_20/obj.C)^(1/4) + ...
             obj.alpha_1 / (2*obj.V_10^(1/2)*obj.A^(1/2)) * (V_1-obj.V_10) - obj.alpha_2 / (4*obj.V_20^(3/4)*obj.C^(1/4)) * (V_2-obj.V_20);
    
    %% Wymuszenia
    F_1in = u(1,:) + obj.F_10;
    F_D = u(2,:) + obj.F_D0;
    
    %% RK4 - Runge-Kutta 4. rzędu
    for i = 2:kk
        if i <= obj.delay + 1
            % RK4 dla nieliniowego układu
            % Krok dla V_1
            k1_V1 = obj.Tp * fun_1(obj.F_10, F_D(i), V_1);
            k2_V1 = obj.Tp * fun_1(obj.F_10, F_D(i), V_1 + 0.5 * k1_V1);
            k3_V1 = obj.Tp * fun_1(obj.F_10, F_D(i), V_1 + 0.5 * k2_V1);
            k4_V1 = obj.Tp * fun_1(obj.F_10, F_D(i), V_1 + k3_V1);
            V_1 = V_1 + (k1_V1 + 2*k2_V1 + 2*k3_V1 + k4_V1) / 6;
    
            % Krok dla V_2
            k1_V2 = obj.Tp * fun_2(V_1, V_2);
            k2_V2 = obj.Tp * fun_2(V_1, V_2 + 0.5 * k1_V2);
            k3_V2 = obj.Tp * fun_2(V_1, V_2 + 0.5 * k2_V2);
            k4_V2 = obj.Tp * fun_2(V_1, V_2 + k3_V2);
            V_2 = V_2 + (k1_V2 + 2*k2_V2 + 2*k3_V2 + k4_V2) / 6;
    
            % RK4 dla układu liniowego
            k1_V1L = obj.Tp * fun_1L(obj.F_10, F_D(i), V_1L);
            k2_V1L = obj.Tp * fun_1L(obj.F_10, F_D(i), V_1L + 0.5 * k1_V1L);
            k3_V1L = obj.Tp * fun_1L(obj.F_10, F_D(i), V_1L + 0.5 * k2_V1L);
            k4_V1L = obj.Tp * fun_1L(obj.F_10, F_D(i), V_1L + k3_V1L);
            V_1L = V_1L + (k1_V1L + 2*k2_V1L + 2*k3_V1L + k4_V1L) / 6;
    
            k1_V2L = obj.Tp * fun_2L(V_1L, V_2L);
            k2_V2L = obj.Tp * fun_2L(V_1L, V_2L + 0.5 * k1_V2L);
            k3_V2L = obj.Tp * fun_2L(V_1L, V_2L + 0.5 * k2_V2L);
            k4_V2L = obj.Tp * fun_2L(V_1L, V_2L + k3_V2L);
            V_2L = V_2L + (k1_V2L + 2*k2_V2L + 2*k3_V2L + k4_V2L) / 6;
        else
            % Uwzględnienie opóźnienia dla F_1in
            delayed_idx = i - obj.delay;
    
            % RK4 dla nieliniowego układu
            k1_V1 = obj.Tp * fun_1(F_1in(delayed_idx), F_D(i), V_1);
            k2_V1 = obj.Tp * fun_1(F_1in(delayed_idx), F_D(i), V_1 + 0.5 * k1_V1);
            k3_V1 = obj.Tp * fun_1(F_1in(delayed_idx), F_D(i), V_1 + 0.5 * k2_V1);
            k4_V1 = obj.Tp * fun_1(F_1in(delayed_idx), F_D(i), V_1 + k3_V1);
            V_1 = V_1 + (k1_V1 + 2*k2_V1 + 2*k3_V1 + k4_V1) / 6;
    
            k1_V2 = obj.Tp * fun_2(V_1, V_2);
            k2_V2 = obj.Tp * fun_2(V_1, V_2 + 0.5 * k1_V2);
            k3_V2 = obj.Tp * fun_2(V_1, V_2 + 0.5 * k2_V2);
            k4_V2 = obj.Tp * fun_2(V_1, V_2 + k3_V2);
            V_2 = V_2 + (k1_V2 + 2*k2_V2 + 2*k3_V2 + k4_V2) / 6;
    
            % RK4 dla układu liniowego
            k1_V1L = obj.Tp * fun_1L(F_1in(delayed_idx), F_D(i), V_1L);
            k2_V1L = obj.Tp * fun_1L(F_1in(delayed_idx), F_D(i), V_1L + 0.5 * k1_V1L);
            k3_V1L = obj.Tp * fun_1L(F_1in(delayed_idx), F_D(i), V_1L + 0.5 * k2_V1L);
            k4_V1L = obj.Tp * fun_1L(F_1in(delayed_idx), F_D(i), V_1L + k3_V1L);
            V_1L = V_1L + (k1_V1L + 2*k2_V1L + 2*k3_V1L + k4_V1L) / 6;
    
            k1_V2L = obj.Tp * fun_2L(V_1L, V_2L);
            k2_V2L = obj.Tp * fun_2L(V_1L, V_2L + 0.5 * k1_V2L);
            k3_V2L = obj.Tp * fun_2L(V_1L, V_2L + 0.5 * k2_V2L);
            k4_V2L = obj.Tp * fun_2L(V_1L, V_2L + k3_V2L);
            V_2L = V_2L + (k1_V2L + 2*k2_V2L + 2*k3_V2L + k4_V2L) / 6;
        end
    
        % Wyznaczanie wysokości słupa cieczy
        h_2e = sqrt(V_2 / obj.C);
        h_2eL = sqrt(obj.V_20 / obj.C) + 1/(2*sqrt(obj.V_20 * obj.C)) * (V_2L - obj.V_20);
    
        % Sprowadzenie wartości do punktu pracy
        y(i) = h_2e - obj.h_20;
        y_L(i) = h_2eL - obj.h_20;
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

    %% Prezentacja wyników
    figure;
    hold on;
    plot(0:obj.Tp:(kk-1)*obj.Tp, y_L, 'b.','LineWidth',2);
    plot(0:obj.Tp:(kk-1)*obj.Tp, y, 'r-','LineWidth',2);
    hold off;
    title('Wysokość słupa cieczy w zbiorniku 2. - h_2(t)');
    legend('h_{2lin}', 'h_2');
    xlabel('t [s]');
    ylabel('h [cm]');
    grid on;
end