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
obiekt = Obiekt('ARX', 'linear');

%% Odpowiedź skokowa
obiekt = obiekt.sopdt();
% obiekt.show_sopdt();

%% Hammerstein
hammerstein = Hammerstein(obiekt.F_10-45, obiekt.F_10+45, obiekt.F_10, ...
                          obiekt.F_D0-20, obiekt.F_D0+20, obiekt.F_D0, ...
                          obiekt.h_20, obiekt.alpha_2, obiekt.n, obiekt.mode, obiekt.s);
hammerstein = hammerstein.fuzzyfication();
% hammerstein.show_fuzzy_system();

%% Wiener
wiener = Wiener(obiekt.a, obiekt.b, obiekt.K, obiekt.n, obiekt.delay, @obiekt.modified_Euler, obiekt.mode, obiekt.s);
wiener = wiener.fuzzyfication();
% wiener.show_fuzzy_system();

%% Wymuszenia
u = [repelem((rand(1, obiekt.kk/250) * 5 - 2.5) * 10, 250); zeros(1, obiekt.kk)];
[y, y_L] = obiekt.modified_Euler(u, obiekt.kk);
% obiekt.show_modified_Euler(u, y, y_L, obiekt.kk);

%% Testowanie
obiekt.diff_eq(u, y)
hammerstein.model(obiekt.a, obiekt.b, u, y, obiekt.K, obiekt.kk, obiekt.delay);
wiener.model(obiekt.a, obiekt.b, u, y, obiekt.K, obiekt.kk, obiekt.delay);