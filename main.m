%% PPMGR
% Autor: Wojciech Rogalski
% Data: 11.03.2024r.

%% Regulacja kaskadowa
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

%% Charakrerystyka statyczna
static_characteristic(F_10, F_D0, h_10, h_20, alpha_1, alpha_2);

%% Wymuszenia
[u] = enforce(kk);

%% Linearyzacja za pomocą zmodyfikowanej metody Eulera
modified_Euler(a, c, alpha_1, alpha_2, F_10, F_D0, h_10, h_20, V_10, V_20, t, kk, tau, Tp, u);

%% Transmitancje G(s) i G(z)
[G_s, G_z] = tf_function(a, c, alpha_1, alpha_2, V_10, V_20, tau, Tp);

%% Przebiegi wartości wyjściowych obliczone przy pomocy G(s), G(z)
tf_continious(G_s, u, t);
tf_discrete(G_z, u, t, kk);

%% Wartość zadana h_2
delta_h = set_value(kk, h_20);

%% Współczynniki do równań różnicowych na h_2
b(1) = [G_z.Numerator{2,1}(2)];
b(2) = [G_z.Numerator{2,1}(3)];

a(1) = [G_z.Denominator{2,1}(2)];
a(2) = [G_z.Denominator{2,1}(3)];

%% PID(z)
% Parametry PID
Kp = 6.07 * 0.45;
Ti = 460/2;
Td = 460/8;

PID_controller(Kp, Ti, Td, a, b, u, delta_h, Tp, kk, tau/Tp, t);

%% Regulator DMC
% Horyzonty
N = 100;
D = N;
Nu = 2;
lambda = 0.2;

% Ograniczenia (x_min = -x_max)
y_max = 10;
u_max = 15;
delta_u_max = 5;

% Odpwiedzi skokowe
s = step(G_z(2,1), Tp*(D-1));

%% Regulacja DMC - analitycznie
DMC_analitic(a, b, N, Nu, D, lambda, s, u, delta_h, y_max, u_max, delta_u_max, kk, tau/Tp, t);

%% Regulacja DMC - numerycznie
DMC_numeric(a, b, N, Nu, D, lambda, s, u, delta_h, y_max, u_max, delta_u_max, kk, tau/Tp, t);

%% Fuzzy static charaacteristic
% Zakres rozmywania y(u)
u_start = 45;
u_end = 135;
[R_U, optimal_params_U] = static_characteristic_y_u(F_D0, alpha_2, u_start, u_end);
% Zakres rozmywania u(y)
y_start = 10;
y_end = 70;
[R_Y, optimal_params_Y] = static_characteristic_u_y(F_D0, alpha_2, y_start, y_end);

%% Fuzzy PID
PID_controller_fuzzy(Kp, Ti, Td, a, b, u, delta_h, Tp, kk, tau/Tp, t, R_U, R_Y, optimal_params_U, optimal_params_Y);

%% Fuzzy DMC - analitic
DMC_analitic_fuzzy(a, b, N, Nu, D, lambda, s, u, delta_h, y_max, u_max, delta_u_max, kk, tau/Tp, t, R_U, R_Y, optimal_params_U, optimal_params_Y);

%% Fuzzy DMC - numeric
DMC_numeric_fuzzy(a, b, N, Nu, D, lambda, s, u, delta_h, y_max, u_max, delta_u_max, kk, tau/Tp, t, R_U, R_Y, optimal_params_U, optimal_params_Y);