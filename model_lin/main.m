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
n = (F_1_end - F_1_start) * 4;

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
% file_name = sprintf('../raport/pictures/static_characteristic.pdf');
% exportgraphics (gcf, file_name);

%% Podział danych statycznych
u = F_1 - F_10;
u_start = F_1_start - F_10;
u_end = F_1_end - F_10;

u_ucz = F_1(1:2:end-1) - F_10;
u_wer = F_1(2:2:end) - F_10;
y_ucz = h_2(1:2:end-1) - h_20;
y_wer = h_2(2:2:end) - h_20;

v_start = -25;
v_end = 25;
v = linspace(v_start, v_end, n);
v_ucz = v(1:2:end-1);
v_wer = v(2:2:end);

%% Fuzzy static charaacteristic - Hammerstein
[R_hammerstein, hammerstein_params] = static_characteristic_y_u(u, u_ucz, u_wer, y_ucz, y_wer, u_start, u_end, n, 'linear');

%% Fuzzy static charaacteristic - Wiener
[R_wiener, wiener_params] = static_characteristic_y_v(v, v_ucz, v_wer, y_ucz, y_wer, v_start, v_end, n, 'linear');

%% Wymuszenia + podział na zbiory danych dynamicznych
u = enforce(kk);
u_ucz = [u(1, 1:kk/2); u(2, 1:kk/2)];
u_wer = [u(1, kk/2+1:end); u(2, kk/2+1:end)];
% file_name = sprintf('../raport/pictures/u_F1.pdf');
% exportgraphics (gcf, file_name);

%% Wyjście OBIEKTU - modified Euler + podział na zbiory danych dynamicznych
% Wyświetlane wartości na wykresach są sprowadzone do pkt. pracy --> wartości przyrostowe
y = modified_Euler(A, C, alpha_1, alpha_2, F_10, F_D0, h_20, V_10, V_20, t, kk, tau, Tp, u);
y_ucz = y(1:kk/2);
y_wer = y(kk/2+1:end);
% file_name = sprintf('../raport/pictures/y_F1.pdf');
% exportgraphics (gcf, file_name);

%% Przebiegi wartości wyjściowych obliczone przy pomocy G(z) - MODEL zlinearyzowany (ARX)
% y = tf_discrete(G_z, u, t, kk);
diff_eq_ARX(a(:,2), b(:,2), u_ucz, u_wer, y_ucz, y_wer, kk, tau/Tp);

%% Przebiegi wartości wyjściowych obliczone przy pomocy G(z) - MODEL zlinearyzowany (OE)
diff_eq_OE(a(:,2), b(:,2), u_ucz, u_wer, y_ucz, y_wer, kk, tau/Tp);

%% Przebiegi wartości wyjściowych obliczone przy pomocy statyki rozmytej - MODEL Hammersteina (ARX)
% Przekazać do funkcji wsp. b podzielone przez dcgain
diff_eq_hammerstein_ARX(a(:,2), b(:,2), dcgain(G_z(2,1)), u_ucz, u_wer, y_ucz, y_wer, hammerstein_params, u_start, u_end, R_hammerstein, kk, tau/Tp, n);

%% Przebiegi wartości wyjściowych obliczone przy pomocy statyki rozmytej - MODEL Hammersteina (OE)
% Przekazać do funkcji wsp. b podzielone przez dcgain
diff_eq_hammerstein_OE(a(:,2), b(:,2), dcgain(G_z(2,1)), u_ucz, u_wer, y_ucz, y_wer, hammerstein_params, u_start, u_end, R_hammerstein, kk, tau/Tp, n);

%% Przebiegi wartości wyjściowych obliczone przy pomocy statyki rozmytej - MODEL Wienera (ARX)
% Przekazać do funkcji wsp. b podzielone przez dcgain
v = diff_eq_wiener_ARX(a(:,2), b(:,2), dcgain(G_z(2,1)), u_ucz, u_wer, y_ucz, y_wer, wiener_params, v_start, v_end, R_wiener, kk, tau/Tp, n);

%% Przebiegi wartości wyjściowych obliczone przy pomocy statyki rozmytej - MODEL Wienera (OE)
% Przekazać do funkcji wsp. b podzielone przez dcgain
diff_eq_wiener_OE(a(:,2), b(:,2), dcgain(G_z(2,1)), u_ucz, u_wer, y_ucz, y_wer, wiener_params, v_start, v_end, R_wiener, kk, tau/Tp, n);