%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kaskadowa regulacja poziomu wody w zbiornikach %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Skrypt porównujący przebiegu h_1(t) i h_2(t) dla modeli: nieliniowego, liniowego, liniowego rozmytego
% Autor: Wojciech Rogalski
% Data: 25.10.2023r.
clear all;

%% Dane
P = 540;
C = 0.85;
alfa_1 = 26;
alfa_2 = 20;
Tp = 50;
steps = 100;

%% Punkt pracy
F_10 = 100;
F_D0 = 30;
h_10 = ((F_10+F_D0)/alfa_1)^2;
h_20 = ((F_10+F_D0)/alfa_2)^2;
V_10 = P * h_10;
V_20 = C * h_20^2;
tau = 100;

%% Alokacja pamięci
F_1 = zeros(1, steps);
F_1in = zeros(1, steps);
F_D = zeros(1,steps);
y_1 = zeros(1, steps);
y_1L = zeros(1, steps);
y_1F = zeros(1, steps);
y_2 = zeros(1, steps);
y_2L = zeros(1, steps);
y_2F = zeros(1, steps);
t = zeros(1, steps);

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

y_1(1) = h_10;
y_1L(1) = h_10;
y_1F(1) = h_10;
y_2(1) = h_20;
y_2L(1) = h_20;
y_2F(1) = h_20;

%% Zbiory rozmyte F_1
V_1F = zeros(1, steps);
V_2F = zeros(1, steps);
h_1F = zeros(1, steps);
h_2F = zeros(1, steps);
V_1F(1:tau/Tp+1) = V_10;
V_2F(1:tau/Tp+1) = V_20;
h_1F(1:tau/Tp+1) = h_10;
h_2F(1:tau/Tp+1) = h_20;

%% Rozmywanie F_1
sets = 5;
n = 1000;
width = 2.5;
F_1Fmax = 120;
F_1Fmin = 80;
F_1F = zeros(1, sets);
F_1F(1) = 80;
F_1F(2) = 90;
F_1F(3) = 100;
F_1F(4) = 110;
F_1F(5) = 120;
x = linspace(F_1Fmin, F_1Fmax, n);
w_arr = cell(1, sets);

for i = 1:sets
    w_arr{i} = gaussmf(x, [width, F_1F(i)]);
end

figure(1);
hold on;
plot(x, w_arr{1});
plot(x, w_arr{2});
plot(x, w_arr{3});
plot(x, w_arr{4});
plot(x, w_arr{5});
hold off;
title('Funkcje przynależności wartości F_{1}');
xlabel('F_{1}');
ylabel('mi(F_{1})');
axis([F_1Fmin F_1Fmax 0 1]);

%% Punkty pracy dla różnych wartości F_1F
V_1F0 = zeros(1, sets);
V_2F0 = zeros(1, sets);
h_1F0 = zeros(1, sets);
h_2F0 = zeros(1, sets);
V_1eLF = zeros(1, sets);
V_2eLF = zeros(1, sets);
h_1eLF = zeros(1, sets);
h_2eLF = zeros(1, sets);
w_1 = zeros(1, n);
w_2 = zeros(1, n);

for i = 1:sets
    h_1F0(i) = ((F_1F(i) + F_D0) / alfa_1)^2;
    h_2F0(i) = ((F_1F(i) + F_D0) / alfa_2)^2;
    V_1F0(i) = P * h_1F0(i);
    V_2F0(i) = C * (h_2F0(i))^2;
end

%% Funkcje różniczkujące
fun_1 = @(h_1, F_1, F_D) F_1 + F_D - alfa_1 * sqrt(h_1);
fun_2 = @(h_1, h_2) alfa_1 * sqrt(h_1) - alfa_2 * sqrt(h_2);
fun_1L = @(h_1, h10, F_1, F_D, F10, FD0) (F_1-F10) + (F_D-FD0) - alfa_1 / (2*sqrt(h10)) * (h_1-h10);
fun_2L = @(h_1, h10, h_2, h20) alfa_1 / (2*sqrt(h10)) * (h_1-h10) - alfa_2 / (2*sqrt(h20)) * (h_2-h20);

F_1(1) = F_10;

%% Dane odpowiadjące wykresom przedstawiającym odpowiedzi układu na zmienę dopływu sterującego F_1
% F_1in(1:25) = F_10+20;
% F_1in(26:50) = F_10-15;
% F_1in(51:steps) = F_10;
% F_D(1:steps) = F_D0;

%% Dane odpowiadjące wykresom przedstawiającym odpowiedzi układu na zmienę dopływu zakłócającego F_1
F_1in(1:steps) = F_10;
F_D(1:25) = F_D0+25;
F_D(26:50) = F_D0-25;
F_D(51:steps) = F_D0;

%% Modified Euler 
for i = 2:steps    
    if(i <= tau/Tp)
        F_1(i) = F_10;
        y_1(i) = h_10;
        y_1L(i) = h_10;
        y_1F(i) = h_10;
        y_2(i) = h_20;
        y_2L(i) = h_20;
        y_2F(i) = h_20;
        t(i) = i*Tp;
        continue;
    else
        F_1(i) = F_1in(i-tau/Tp);
    end
    
    %% Część nieliniowa
    V_1 = V_1e + Tp * fun_1(h_1e, F_1(i-1), F_D(i-1));
    V_2 = V_2e + Tp * fun_2(h_1e, h_2e);
    h_1 = V_1 / P;
    h_2 = sqrt(V_2 / C);

    V_1e = V_1e + 1/2 * Tp * (fun_1(h_1e, F_1(i-1), F_D(i-1)) + fun_1(h_1, F_1(i), F_D(i)));  
    V_2e = V_2e + 1/2 * Tp * (fun_2(h_1e, h_2e) + fun_2(h_1, h_2));
    h_1e = V_1e / P;
    h_2e = sqrt(V_2e / C);

    %% Część liniowa
    V_1L = V_1eL + Tp * fun_1L(h_1eL, h_10, F_1(i-1), F_D(i-1), F_10, F_D0);
    V_2L = V_2eL + Tp * fun_2L(h_1eL, h_10, h_2eL, h_20);
    h_1L= V_1L / P;
    h_2L = sqrt(V_2L / C);

    V_1eL = V_1eL + 1/2 * Tp * (fun_1L(h_1eL, h_10, F_1(i-1), F_D(i-1), F_10, F_D0) + fun_1L(h_1L, h_10, F_1(i), F_D(i), F_10, F_D0));
    V_2eL = V_2eL + 1/2 * Tp * (fun_2L(h_1eL, h_10, h_2eL, h_20) + fun_2L(h_1L, h_10, h_2L, h_20));
    h_1eL = V_1eL / P;
    h_2eL = sqrt(V_2eL / C);

    %% Część rozmyta
    number_1 = round((F_1(i) - F_1Fmin) * (n - 0)/(F_1Fmax - F_1Fmin) + 0);
    number_2 = round((F_1(i-1) - F_1Fmin) * (n - 0)/(F_1Fmax - F_1Fmin) + 0);
    for j = 1:sets
        w_1(j) = w_arr{j}(number_1);
        w_2(j) = w_arr{j}(number_2);

        V_1LF = V_1F(i-1) + Tp * fun_1L(h_1F(i-1), h_1F0(j), F_1(i-1), F_D(i-1), F_1F(j), F_D0);
        V_2LF = V_2F(i-1) + Tp * fun_2L(h_1F(i-1), h_1F0(j), h_2F(i-1), h_2F0(j));
        h_1LF = V_1LF / P;
        h_2LF = sqrt(V_2LF / C);
    
        V_1eLF(j) = V_1F(i-1) + 1/2 * Tp * (fun_1L(h_1F(i-1), h_1F0(j), F_1(i-1), F_D(i-1), F_1F(j), F_D0) + fun_1L(h_1LF, h_1F0(j), F_1F(j), F_D(i), F_1F(j), F_D0));
        V_2eLF(j) = V_2F(i-1) + 1/2 * Tp * (fun_2L(h_1F(i-1), h_1F0(j), h_2F(i-1), h_2F0(j)) + fun_2L(h_1LF, h_1F0(j), h_2LF, h_2F0(j)));
        h_1eLF(j) = V_1eLF(j) / P;
        h_2eLF(j) = sqrt(V_2eLF(j) / C);
    end

        V_1_sum = 0;
        V_2_sum = 0;
        h_1_sum = 0;
        h_2_sum = 0;
        w_sum_1 = 0;
        w_sum_2 = 0;

    for j = 1:sets
        V_1_sum = V_1_sum + V_1eLF(j) * w_1(j);
        V_2_sum = V_2_sum + V_2eLF(j) * w_2(j);
        h_1_sum = h_1_sum + h_1eLF(j) * w_1(j);
        h_2_sum = h_2_sum + h_2eLF(j) * w_2(j);
        w_sum_1 = w_sum_1 + w_1(j);
        w_sum_2 = w_sum_2 + w_2(j);
    end

    V_1F(i) = V_1_sum / w_sum_1;
    V_2F(i) = V_2_sum / w_sum_2;
    h_1F(i) = h_1_sum / w_sum_1;
    h_2F(i) = h_2_sum / w_sum_2;

    y_1(i) = h_1e;
    y_1L(i) = h_1eL;
    y_1F(i) = h_1F(i);
    y_2(i) = h_2e;
    y_2L(i) = h_2eL;
    y_2F(i) = h_2F(i);
    t(i) = i*Tp;
end

%% Prezentacja wyników
figure(2);
hold on;
plot(t, y_1, 'r--','LineWidth',2);
plot(t, y_1L, 'b','LineWidth',2);
plot(t, y_1F, 'g-.', 'LineWidth',2);
hold off;
title('h_{1}(t)');
legend('h_1', 'h_{1lin}', 'h_{1fuzzy}');
xlabel('t [s]');
ylabel('h [cm]');
grid on;

figure(3);
hold on;
plot(t, y_2, 'g--','LineWidth',2);
plot(t, y_2L, 'y','LineWidth',2);
plot(t, y_2F, 'm-.', 'LineWidth',2);
hold off;
title('h_{2}(t)');
legend('h_2', 'h_{2lin}', 'h_{2fuzzy}');
xlabel('t [s]');
ylabel('h [cm]');
grid on;