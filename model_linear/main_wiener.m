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
% n = (F_1_end - F_1_start) * 4;
n = 5000;

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

%% Podział danych statycznych
v_start = -25;
v_end = 25;
u = zeros(2, n);
u(1, :) = linspace(-45, 45, n);
u(2, :) = 0;

%% Wyjście OBIEKTU - modified Euler + podział na zbiory danych dynamicznych
% Wyświetlane wartości na wykresach są sprowadzone do pkt. pracy --> wartości przyrostowe
y = modified_Euler(A, C, alpha_1, alpha_2, F_10, F_D0, h_20, V_10, V_20, t, n, tau, Tp, u);
y_ucz_stat = y(500:2:end-1);
y_wer_stat = y(501:2:end);
% file_name = sprintf('../raport/pictures/y_F1.pdf');
% exportgraphics (gcf, file_name);

%% Fuzzy static charaacteristic - Wiener
v = v_ARX(a(:,2), b(:,2), dcgain(G_z(2,1)), u, y, n, tau/Tp);
v_ucz_stat = v(500:2:end-1);
v_wer_stat = v(501:2:end);

%% Wiener params
[R_wiener, wiener_params] = static_characteristic_y_v(v(500:end), v_ucz_stat, v_wer_stat, y_ucz_stat, y_wer_stat, v_start, v_end, n-500, 'linear');

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

%% Przebiegi wartości wyjściowych obliczone przy pomocy statyki rozmytej - MODEL Wienera (ARX)
% Przekazać do funkcji wsp. b podzielone przez dcgain
diff_eq_wiener_ARX(a(:,2), b(:,2), dcgain(G_z(2,1)), u_ucz, u_wer, y_ucz, y_wer, wiener_params, v_start, v_end, R_wiener, kk, tau/Tp, n-500);

%% Podział danych statycznych
v_start = -30;
v_end = 30;
u = zeros(2,n);
u(1,:) = linspace(-45, 45, n);
u(2,:) = 0;

%% Wyjście OBIEKTU - modified Euler + podział na zbiory danych dynamicznych
% Wyświetlane wartości na wykresach są sprowadzone do pkt. pracy --> wartości przyrostowe
y = modified_Euler(A, C, alpha_1, alpha_2, F_10, F_D0, h_20, V_10, V_20, t, n, tau, Tp, u);
y_ucz_stat = y(500:2:end-1);
y_wer_stat = y(501:2:end);
% file_name = sprintf('../raport/pictures/y_F1.pdf');
% exportgraphics (gcf, file_name);

%% Fuzzy static charaacteristic - Wiener
v = v_OE(a(:,2), b(:,2), dcgain(G_z(2,1)), u, y, n, tau/Tp);
v_ucz_stat = v(500:2:end-1);
v_wer_stat = v(501:2:end);

%% Wiener params
[R_wiener, wiener_params] = static_characteristic_y_v(v(500:end), v_ucz_stat, v_wer_stat, y_ucz_stat, y_wer_stat, v_start, v_end, n-500, 'linear');

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

%% Przebiegi wartości wyjściowych obliczone przy pomocy statyki rozmytej - MODEL Wienera (OE)
% Przekazać do funkcji wsp. b podzielone przez dcgain
diff_eq_wiener_OE(a(:,2), b(:,2), dcgain(G_z(2,1)), u_ucz, u_wer, y_ucz, y_wer, wiener_params, v_start, v_end, R_wiener, kk, tau/Tp, n-500);

%% Przebiegi wartości wyjściowych obliczone przy pomocy G(z) - MODEL zlinearyzowany (ARX)
% y = tf_discrete(G_z, u, t, kk);
diff_eq_ARX(a(:,2), b(:,2), u_ucz, u_wer, y_ucz, y_wer, kk, tau/Tp);

%% Przebiegi wartości wyjściowych obliczone przy pomocy G(z) - MODEL zlinearyzowany (OE)
diff_eq_OE(a(:,2), b(:,2), u_ucz, u_wer, y_ucz, y_wer, kk, tau/Tp);