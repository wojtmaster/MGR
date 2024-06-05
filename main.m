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
A = 540;
C = 0.85;
alpha_1 = 26;
alpha_2 = 20;

% Punkt pracy
F_10 = 90;
F_D0 = 30;
h_10 = ((F_10+F_D0)/alpha_1)^2;
h_20 = ((F_10+F_D0)/alpha_2)^2;
V_10 = A * h_10;
V_20 = C * h_20^2;

% Opóźnienie
tau = 100;
% Okres próbkowania
Tp = 20;
% Próbki dyskretne
kk = 5000;
% Wektor czasu
t = 0:Tp:(kk-1)*Tp;
% Zakres wyznaczanej charakterystyki statycznej
F_1_start = 45;
F_1_end = 135;
% Liczba podziałów wartości sterowania F_1 (do charakterystyki statycznej)
n = (F_1_end - F_1_start) * 2;

%% Transmitancje G(z)
[G_z] = tf_function(A, C, alpha_1, alpha_2, V_10, V_20, tau, Tp);
% Współczynniki do równań różnicowych na h_1
a(1,1) = [G_z.Denominator{1,1}(2)];
b(tau/Tp+1,1) = [G_z.Numerator{1,1}(2)];
% Współczynniki do równań różnicowych na h_2
b(tau/Tp+1,2) = [G_z.Numerator{2,1}(2)];
b(tau/Tp+2,2) = [G_z.Numerator{2,1}(3)];
a(1,2) = [G_z.Denominator{2,1}(2)];
a(2,2) = [G_z.Denominator{2,1}(3)];

%% Charakrerystyka statyczna
[h_2, F_1] = static_characteristic(F_1_start, F_1_end, F_D0, alpha_2, n);

%% Podział danych statycznych
F_1_ucz = F_1(1:2:end-1);
F_1_wer = F_1(2:2:end);
h_2_ucz = h_2(1:2:end-1);
h_2_wer = h_2(2:2:end);

%% Fuzzy static charaacteristic - linear
[R, optimal_params] = static_characteristic_y_u(F_1, F_1_ucz, F_1_wer, h_2_ucz, h_2_wer, F_1_start, F_1_end, n, 'linear');

%% Fuzzy static charaacteristic - nonlinear
% [R, optimal_params] = static_characteristic_y_u(F_1, F_1_ucz, F_1_wer, h_2_ucz, h_2_wer, F_1_start, F_1_end, n, 'nonlinear');

%% Wymuszenia + podział na zbiory danych dynamicznych
u = enforce(kk);
u_ucz = [u(1, 1:kk/2); u(2, 1:kk/2)];
u_wer = [u(1, kk/2+1:end); u(2, kk/2+1:end)];

%% Wyjście OBIEKTU - modified Euler + podział na zbiory danych dynamicznych
% Wyświetlane wartości na wykresach są sprowadzone do pkt. pracy --> wartości przyrostowe
y = modified_Euler(A, C, alpha_1, alpha_2, F_10, F_D0, h_20, V_10, V_20, t, kk, tau, Tp, u);
y_ucz = y(1:kk/2);
y_wer = y(kk/2+1:end);

%% Przebiegi wartości wyjściowych obliczone przy pomocy G(z) - MODEL zlinearyzowany (ARX)
% y = tf_discrete(G_z, u, t, kk);
diff_eq_ARX(a(:,2), b(:,2), u_ucz, u_wer, y_ucz, y_wer, kk, tau/Tp);

%% Przebiegi wartości wyjściowych obliczone przy pomocy G(z) - MODEL zlinearyzowany (OE)
diff_eq_OE(a(:,2), b(:,2), u_ucz, u_wer, y_ucz, y_wer, kk, tau/Tp);

%% Przebiegi wartości wyjściowych obliczone przy pomocy statyki rozmytej - MODEL Hammersteina (ARX)
% Przekazać do funkcji wsp. b podzielone przez dcgain
diff_eq_hammerstein_ARX(a(:,2), b(:,2), dcgain(G_z(2,1)), u_ucz, u_wer, y_ucz, y_wer, optimal_params, F_10, h_20, F_1_start, F_1_end, R, kk, tau/Tp, n);

%% Przebiegi wartości wyjściowych obliczone przy pomocy statyki rozmytej - MODEL Hammersteina (OE)
% Przekazać do funkcji wsp. b podzielone przez dcgain
diff_eq_hammerstein_OE(a(:,2), b(:,2), dcgain(G_z(2,1)), u_ucz, u_wer, y_ucz, y_wer, optimal_params, F_10, h_20, F_1_start, F_1_end, R, kk, tau/Tp, n);

%% Przebiegi wartości wyjściowych obliczone przy pomocy statyki rozmytej - MODEL Wienera (ARX)
% Przekazać do funkcji wsp. b podzielone przez dcgain
v = diff_eq_wiener_ARX(a(:,2), b(:,2), dcgain(G_z(2,1)), u_ucz, u_wer, y_ucz, y_wer, optimal_params, F_10, h_20, F_1_start, F_1_end, R, kk, tau/Tp, n);

%% Przebiegi wartości wyjściowych obliczone przy pomocy statyki rozmytej - MODEL Wienera (OE)
% Przekazać do funkcji wsp. b podzielone przez dcgain
diff_eq_wiener_OE(a(:,2), b(:,2), dcgain(G_z(2,1)), u_ucz, u_wer, y_ucz, y_wer, optimal_params, F_10, h_20, F_1_start, F_1_end, R, kk, tau/Tp, n);

%% DMC-SL
% Horyzonty
D = 100;
N = 25;
Nu = 3;
lambda = 0.25;
% Ograniczenia (x_min = -x_max)
y_max = 50;
u_max = 150;
delta_u_max = 50;
y_zad = set_value(kk, h_20);
DMC_SL(a, b, N, Nu, D, lambda, u, y_zad, y_max, u_max, delta_u_max, A, C, alpha_1, alpha_2, h_10, h_20, kk, tau/Tp, Tp)
s = step(G_z(2,1), Tp*(D-1));
DMC_analitic(a(:,2), b(:,2), N, Nu, D, lambda, s, u, y_zad, y_max, u_max, delta_u_max, kk, tau/Tp);