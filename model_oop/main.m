%% MGR
% Autor: Wojciech Rogalski
% Data: 15.12.2024r.
% Tytuł: Porównanie modeli Hammersteina i Winera w regulacji rozmytej
clear all;
set(groot, 'DefaultTextFontName', 'Arial');
set(groot, 'DefaultAxesFontName', 'Arial');
set(groot, 'DefaultTextFontSize', 12);
set(groot, 'DefaultAxesFontSize', 10);

%% Dane
obiekt = Obiekt('OE', 'nonlinear');

%% Odpowiedź skokowa na wymuszenie F_1
u = [ones(1, obiekt.kk); zeros(1, obiekt.kk)];
[y, ~] = obiekt.rk4(u, obiekt.kk);
[a, b, K, G_z] = obiekt.sopdt(0.6025, obiekt.tau, 200, 10);
% obiekt.show_sopdt(y, G_z);

%% Odpowiedź skokowa na zakłócenie F_D
u = [zeros(1, obiekt.kk); ones(1, obiekt.kk)];
[y, ~] = obiekt.rk4(u, obiekt.kk);
[a_disturbance, b_disturbance, K_disturbance, G_z_disturbance] = obiekt.sopdt(0.6025, 0, 200, 10);
% obiekt.show_sopdt(y, G_z_disturbance);

% Z racji na to, że 
% a = a_disturbance
% b = b_disturbance 
% K = K_disturbance
% dalej przyjęto jednolitą nomenklaturę a, b, K

%% Hammerstein
hammerstein = Hammerstein(obiekt.F_10-60, obiekt.F_10+60, obiekt.F_10, ...
                          obiekt.F_D0-20, obiekt.F_D0+20, obiekt.F_D0, ...
                          obiekt.h_20, obiekt.alpha_2, obiekt.n, obiekt.mode, obiekt.s);
hammerstein = hammerstein.fuzzyfication();
hammerstein = hammerstein.fuzzyfication_disturbance();

% hammerstein.show_fuzzy_system(hammerstein.u, hammerstein.u_ucz, hammerstein.u_wer, hammerstein.y_ucz, hammerstein.y_mod_ucz, hammerstein.y_wer, hammerstein.y_mod_wer, hammerstein.R, 'F_{1}');
% hammerstein.show_fuzzy_system(hammerstein.u_disturbance, hammerstein.u_ucz_disturbance, hammerstein.u_wer_disturbance, hammerstein.y_ucz_disturbance, hammerstein.y_mod_ucz_disturbance, hammerstein.y_wer_disturbance, hammerstein.y_mod_wer_disturbance, hammerstein.R_disturbance, 'F_{D}');
% hammerstein.show_static_characteristic(hammerstein.u, hammerstein.y, 'F_{1}');
% hammerstein.show_static_characteristic(hammerstein.u_disturbance, hammerstein.y_disturbance, 'F_{D}');

%% Wiener
wiener = Wiener(a, b, K, obiekt.n, obiekt.delay, @obiekt.rk4, obiekt.mode, obiekt.s);
wiener = wiener.fuzzyfication();
wiener = wiener.fuzzyfication_disturbance();
% wiener.show_fuzzy_system(wiener.v, wiener.v_ucz, wiener.v_wer, wiener.y_ucz, wiener.y_mod_ucz, wiener.y_wer, wiener.y_mod_wer, wiener.R, 'F_{1}');
% wiener.show_fuzzy_system(wiener.v_disturbance, wiener.v_ucz_disturbance, wiener.v_wer_disturbance, wiener.y_ucz_disturbance, wiener.y_mod_ucz_disturbance, wiener.y_wer_disturbance, wiener.y_mod_wer_disturbance, wiener.R_disturbance, 'F_{D}');
% wiener.show_static_characteristic(wiener.v, wiener.y, 'F_{1}');
% wiener.show_static_characteristic(wiener.v_disturbance, wiener.y_disturbance, 'F_{D}');

%% Wymuszenia
u = [repelem((rand(1, obiekt.kk/500) * 100 - 50), 500); repelem((rand(1, obiekt.kk/1000) * 20 - 10), 1000)];
[y, y_L] = obiekt.rk4(u, obiekt.kk);
obiekt.show_rk4(u, y, y_L, obiekt.kk);

%% Testowanie
obiekt.diff_eq(a, b, u, y)
hammerstein.model(a, b, K, u, y, obiekt.kk, obiekt.delay);
wiener.model(a, b, K, u, y, obiekt.kk, obiekt.delay);