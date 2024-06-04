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

y(1) = 0;
y_L(1) = 0;

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
    if i <= tau/Tp     
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
plot(t, y, 'r.','LineWidth',2);
hold off;
title('Wysokość słupa cieczy w zbiorniku 2. - h_2(t)');
legend('h_{2lin}', 'h_2');
xlabel('t [s]');
ylabel('h [cm]');
grid on;
end