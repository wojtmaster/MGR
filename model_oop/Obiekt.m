classdef Obiekt
    properties

    end

    methods
        function obj = Obiekt(obj)

        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [u] = enforce(kk)
    u = zeros(2,kk);
    % u(1, :) = F_1
    for i = 1:kk/250
        u(1, (i-1)*250+1 : i*250) = (rand()*5 - 2.5) * 10;
    end
    for i = 1:kk
        if(u(1,i) < -90)
            u(1,i) = -90;
        end
    end
    
    % u(2, :) = F_D
    u(2, :) = 0;
    % u(2, 1 : 0.25*kk) = 0.2;
    % u(2, 0.25*kk+1 : 0.5*kk) = -0.2;
    % u(2, 0.5*kk+1 : 0.75*kk) = -0.2;
    % u(2, 0.75*kk+1 : end) = 0.1;
    % for i = 1:kk
    %     if(u(2,i) < -30)
    %         u(2,i) = -30;
    %     end
    % end

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [y] = modified_Euler(a, c, alpha_1, alpha_2, F_10, F_D0, h_20, V_10, V_20, t, kk, tau, Tp, u)

%% Skrypt przedstawiający zmodyfikowaną metodę Eulera

%% Alokacja pamięci
y = zeros(1, kk);
y_L = zeros(1, kk);

%% Warunki początkowe
V_1 = V_10;
V_2 = V_20;
V_1e = V_10;
V_2e = V_20;
V_1L = V_10;
V_2L = V_20;
V_1eL = V_10;
V_2eL = V_20;

%% Funkcje
fun_1 = @(F_1, F_D, V_1) F_1 + F_D - alpha_1 * (V_1/a)^(1/2);
fun_1L = @(F_1, F_D, V_1) F_1 + F_D - alpha_1 * (V_10/a)^(1/2) - alpha_1 / (2*V_10^(1/2)*a^(1/2)) * (V_1-V_10);
fun_2 = @(V_1, V_2) alpha_1 * (V_1/a)^(1/2) - alpha_2 * (V_2/c)^(1/4); 
fun_2L = @(V_1, V_2) alpha_1 * (V_10/a)^(1/2) - alpha_2 * (V_20/c)^(1/4) + ...
         alpha_1 / (2*V_10^(1/2)*a^(1/2)) * (V_1-V_10) - alpha_2 / (4*V_20^(3/4)*c^(1/4)) * (V_2-V_20);

%% Wymuszenia
F_1in = u(1,:) + F_10;
F_D = u(2,:) + F_D0;

%% Modified Euler 
for i = 2:kk
    if i <= tau/Tp+1    
        V_1 = V_1 + Tp * fun_1(F_10, F_D(i), V_1);
        V_2 = V_2 + Tp * fun_2(V_1, V_2);

        V_1e = V_1e + 1/2 * Tp * (fun_1(F_10, F_D(i), V_1e) + fun_1(F_10, F_D(i), V_1));
        V_2e = V_2e + 1/2 * Tp * (fun_2(V_1e, V_2e) + fun_2(V_1, V_2));

        V_1L = V_1L + Tp * fun_1L(F_10, F_D(i), V_1L);
        V_2L = V_2L + Tp * fun_2L(V_1L, V_2L);

        V_1eL = V_1eL + 1/2 * Tp * (fun_1L(F_10, F_D(i), V_1eL) + fun_1L(F_10, F_D(i), V_1L));
        V_2eL = V_2eL + 1/2 * Tp * (fun_2L(V_1eL, V_2eL) + fun_2L(V_1L, V_2L));
    else                                      
        V_1 = V_1 + Tp * fun_1(F_1in(i - tau/Tp), F_D(i), V_1);
        V_2 = V_2 + Tp * fun_2(V_1, V_2);

        V_1e = V_1e + 1/2 * Tp * (fun_1(F_1in(i - tau/Tp), F_D(i), V_1e) + fun_1(F_1in(i - tau/Tp), F_D(i), V_1));
        V_2e = V_2e + 1/2 * Tp * (fun_2(V_1e, V_2e) + fun_2(V_1, V_2));

        V_1L = V_1L + Tp * fun_1L(F_1in(i - tau/Tp), F_D(i), V_1L);
        V_2L = V_2L + Tp * fun_2L(V_1L, V_2L);

        V_1eL = V_1eL + 1/2 * Tp * (fun_1L(F_1in(i - tau/Tp), F_D(i), V_1eL) + fun_1L(F_1in(i - tau/Tp), F_D(i), V_1L));
        V_2eL = V_2eL + 1/2 * Tp * (fun_2L(V_1eL, V_2eL) + fun_2L(V_1L, V_2L));
    end
        h_2e = sqrt(V_2e / c);
        h_2eL = sqrt(V_20 / c) + 1/(2*sqrt(V_20 * c)) * (V_2eL-V_20);
        
        % Sprowadzenie wartości do punktu pracy
        y(i) = h_2e - h_20;
        y_L(i) = h_2eL - h_20;
end

%% Prezentacja wyników
figure;
hold on;
plot(t, y_L, 'b.','LineWidth',2);
plot(t, y, 'r-','LineWidth',2);
hold off;
title('Wysokość słupa cieczy w zbiorniku 2. - h_2(t)');
legend('h_{2lin}', 'h_2');
xlabel('t [s]');
ylabel('h [cm]');
grid on;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [a, b, K] = sopdt(y, t, tau, Tp)
    T_0 = tau;
    K_0 = 0.6025;
    T_1 = 212;
    T_2 = 15;
    num = K_0;
    den =conv([T_1 1], [T_2 1]);
    
    G_s = tf(num, den, 'InputDelay', T_0);
    
    G_z = c2d(G_s, Tp, 'zoh');
    G_z.Variable = 'z^-1';
    
    a(1:2) = G_z.Denominator{1}(2:end);
    b(1:2) = G_z.Numerator{1}(2:end);
    K = dcgain(G_z);

    %% Prezentacja wyników SOPDT
    figure;
    plot(t, y, 'r-','LineWidth',2);
    hold on;
    step(G_z);
    grid on;
    legend('h', 'h^{mod}', 'Location', 'northwest');
    % file_name = sprintf('../raport/pictures/model_sopdt.pdf');
    % exportgraphics (gcf, file_name);
end
