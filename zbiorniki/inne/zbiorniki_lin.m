%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Regulacja przepływu wody między dwoma zbiornikami %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Porównanie odpowiedzi obiektu liniowego i nieliniowego
clear all;

%% Dane
P = 540;
C = 0.85;
alfa_1 = 26;
alfa_2 = 20;
Tp = 10;
steps = 100;

%% Punkt pracy
F_10 = 90;
F_D0 = 30;
h_10 = ((F_10+F_D0)/alfa_1)^2;
h_20 = ((F_10+F_D0)/alfa_2)^2;
V_10 = P * h_10;
V_20 = C * h_20^2;
tau = 100;

%% Alokacja pamięci
F_1 = zeros(1, steps);
F_D = zeros(1,steps);
t = zeros(1, steps);
y_1 = zeros(1, steps);
y_1L = zeros(1, steps);
y_2 = zeros(1, steps);
y_2L = zeros(1, steps);

%% Warunki początkowe
V_1 = V_10;
V_2 = V_20;
h_1 = h_10;
h_2 = h_20;
V_1e = V_10;
V_2e = V_20;
h_1e = h_10;
h_2e = h_20;
V_1L = V_10;
V_2L = V_20;
h_1L = h_10;
h_2L = h_20;
V_1eL = V_10;
V_2eL = V_20;
h_1eL = h_10;
h_2eL = h_20;

F_1(1) = F_10;
F_D(1) = F_D0;
y_1(1) = h_10;
y_1L(1) = h_10;
y_2(1) = h_20;
y_2L(1) = h_20;

%% Modified Euler 
for i = 2:steps
    if i <= 1/Tp*tau     
        F_1(i) = 95;    
        F_D(i) = F_D0;

        V_1 = V_1 + Tp * fun_1(h_1, F_10, F_D(i), alfa_1);
        V_2 = V_2 + Tp * fun_2(h_1, h_2, alfa_1, alfa_2);
        h_1 = V_1 / P;
        h_2 = sqrt(V_2 / C);

        V_1e = V_1e + 1/2 * Tp * (fun_1(h_1e, F_10, F_D(i-1), alfa_1) + fun_1(h_1, F_10, F_D(i), alfa_1));
        V_2e = V_2e + 1/2 * Tp * (fun_2(h_1e, h_2e, alfa_1, alfa_2) + fun_2(h_1, h_2, alfa_1, alfa_2));
        h_1e = V_1e / P;
        h_2e = sqrt(V_2e / C);

        V_1L = V_1L + Tp * fun_1L(h_1L, F_10, F_D(i-1), alfa_1, h_10);
        V_2L = V_2L+ Tp * fun_2L(h_1L, h_2L, alfa_1, alfa_2, h_10, h_20);
        h_1L= V_1L / P;
        h_2L = sqrt(V_2L / C);

        V_1eL = V_1eL + 1/2 * Tp * (fun_1L(h_1eL, F_10, F_D(i-1), alfa_1, h_10) + fun_1L(h_1L, F_10, F_D(i), alfa_1, h_10));
        V_2eL = V_2eL + 1/2 * Tp * (fun_2L(h_1eL, h_2eL, alfa_1, alfa_2, h_10, h_20) + fun_2L(h_1L, h_2L, alfa_1, alfa_2, h_10, h_20));
        h_1eL = V_1eL / P;
        h_2eL = sqrt(V_2eL / C);
        
        y_1(i) = h_1e;
        y_1L(i) = h_1eL;
        y_2(i) = h_2e;
        y_2L(i) = h_2eL;
        t(i) = (i-1)*Tp;
    else                                   
        F_1(i) = 85;        
        F_D(i) = F_D0;      

        V_1 = V_1 + Tp * fun_1(h_1, F_1(i - 1/Tp*tau), F_D(i-1), alfa_1);
        V_2 = V_2 + Tp * fun_2(h_1, h_2, alfa_1, alfa_2);
        h_1 = V_1 / P;
        h_2 = sqrt(V_2 / C);

        V_1e = V_1e + 1/2 * Tp * (fun_1(h_1e, F_1(i - 1/Tp*tau), F_D(i-1), alfa_1) + fun_1(h_1, F_1(i - 1/Tp*tau + 1), F_D(i), alfa_1));
        V_2e = V_2e + 1/2 * Tp * (fun_2(h_1e, h_2e, alfa_1, alfa_2) + fun_2(h_1, h_2, alfa_1, alfa_2));
        h_1e = V_1e / P;
        h_2e = sqrt(V_2e / C);

        V_1L = V_1L + Tp * fun_1L(h_1L, F_1(i - 1/Tp*tau), F_D(i-1), alfa_1, h_10);
        V_2L = V_2L + Tp * fun_2L(h_1L, h_2L, alfa_1, alfa_2, h_10, h_20);
        h_1L = V_1L / P;
        h_2L = sqrt(V_2L / C);

        V_1eL = V_1eL + 1/2 * Tp * (fun_1L(h_1eL, F_1(i - 1/Tp*tau), F_D(i-1), alfa_1, h_10) + fun_1L(h_1L, F_1(i - 1/Tp*tau + 1), F_D(i), alfa_1, h_10));
        V_2eL = V_2eL + 1/2 * Tp * (fun_2L(h_1eL, h_2eL, alfa_1, alfa_2, h_10, h_20) + fun_2L(h_1L, h_2L, alfa_1, alfa_2, h_10, h_20));
        h_1eL = V_1eL / P;
        h_2eL = sqrt(V_2eL / C);

        y_1(i) = h_1e;
        y_1L(i) = h_1eL;
        y_2(i) = h_2e;
        y_2L(i) = h_2eL;
        t(i) = (i-1)*Tp;
    end
end

%% Rysowanie wykresu
subplot(2,1,1);
hold on;
plot(t, y_1, 'r','LineWidth',2);
plot(t, y_1L, 'b--','LineWidth',2);
hold off;
title('h_{1}(t)');
legend('h_1', 'h_{1lin}');
xlabel('t [s]');
ylabel('h [cm]');
grid on;

subplot(2,1,2);
hold on;
plot(t, y_2, 'g','LineWidth',2);
plot(t, y_2L, 'y--','LineWidth',2);
hold off;
title('h_{2}(t)');
legend('h_2', 'h_{2lin}');
xlabel('t [s]');
ylabel('h [cm]');
grid on;