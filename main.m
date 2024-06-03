%% PPMGR
% Autor: Wojciech Rogalski
% Data: 11.03.2024r.
% Tytuł: Regulacja kaskadowa
clear all;
set(groot, 'DefaultTextFontName', 'Arial');
set(groot, 'DefaultAxesFontName', 'Arial');
set(groot, 'DefaultTextFontSize', 12);
set(groot, 'DefaultAxesFontSize', 10);

%% Dane
a = 540;
c = 0.85;
alpha_1 = 26;
alpha_2 = 20;

% Punkt pracy
F_10 = 90;
F_D0 = 30;
h_10 = ((F_10+F_D0)/alpha_1)^2;
h_20 = ((F_10+F_D0)/alpha_2)^2;
V_10 = a * h_10;
V_20 = c * h_20^2;

% Opóźnienie
tau = 100;
% Okres próbkowania
Tp = 20;
% Próbki dyskretne
kk = 500;
% Wektor czasu
t = 0:Tp:(kk-1)*Tp;
% Zakres wyznaczanej charakterystyki statycznej
F_1_start = 45;
F_1_end = 135;
% Liczba podziałów wartości sterowania F_1 (do charakterystyki statycznej)
n = 200;

%% Transmitancje G(z)
[G_z] = tf_function(a, c, alpha_1, alpha_2, V_10, V_20, tau, Tp);

%% Charakrerystyka statyczna
[h_2, F_1] = static_characteristic(F_1_start, F_1_end, F_D0, alpha_2, n);

%% Fuzzy static charaacteristic - linear
[R, optimal_params] = static_characteristic_y_u(F_1, h_2, n, 'linear');

%% Fuzzy static charaacteristic - nonlinear
% [R, optimal_params] = static_characteristic_y_u(F_D0, alpha_2, F_1_start, F_1_end, n, 'nonlinear');

%% Wymuszenia
[u] = enforce(kk);

%% Linearyzacja za pomocą zmodyfikowanej metody Eulera
modified_Euler(a, c, alpha_1, alpha_2, F_10, F_D0, h_10, h_20, V_10, V_20, t, kk, tau, Tp, u);

%% Przebiegi wartości wyjściowych obliczone przy pomocy G(s), G(z)
tf_discrete(G_z, u, t, kk);

%% Wartość zadana h_2
y_zad = set_value(kk, h_20);

%% Współczynniki do równań różnicowych na h_2
b(1) = [G_z.Numerator{1}(2)];
b(2) = [G_z.Numerator{1}(3)];

a(1) = [G_z.Denominator{1}(2)];
a(2) = [G_z.Denominator{1}(3)];

%% Regulator DMC
% Horyzonty
D = 100;
N = 25;
Nu = 3;
lambda = 0.25;

% Ograniczenia (x_min = -x_max)
y_max = 10;
u_max = 15;
delta_u_max = 5;

% Odpwiedzi skokowe
s = step(G_z(1), Tp*(D-1));

%% Regulacja DMC - analitycznie
DMC_analitic(a, b, N, Nu, D, lambda, s, u, y_zad, y_max, u_max, delta_u_max, kk, tau/Tp);

%% Regulacja DMC - numerycznie
DMC_numeric(a, b, N, Nu, D, lambda, s, u, y_zad, y_max, u_max, delta_u_max, kk, tau/Tp);

%% Regulacja DMC-SL - analitycznie




%% Fuzzy DMC - analitic
DMC_analitic_fuzzy(F_10, h_20, a, b, dcgain(G_z(2,1)), N, Nu, D, lambda, s, u, y_zad, y_max, u_max, delta_u_max, kk, tau/Tp, R, optimal_params, F_1_start, F_1_end, n);

%% Fuzzy DMC - numeric
% DMC_numeric_fuzzy(F_10, h_20, a, b, N, Nu, D, lambda, s, u, y_zad, y_max, u_max, delta_u_max, kk, tau/Tp, R, optimal_params, F_1_start, F_1_end, n);