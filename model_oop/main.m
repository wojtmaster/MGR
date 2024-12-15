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
% Podział danych
n = 500;
% Mode - ARX / OE
mode = 'OE';
% Następniki - linear / nonlinear
s = 'nonlinear';

%% Odpowiedź skokowa
u = [ones(1, kk); zeros(1, kk)];
y = modified_Euler(A, C, alpha_1, alpha_2, F_10, F_D0, h_20, V_10, V_20, t, kk, tau, Tp, u);

%% Model inercjalny SOPDT na podstawie odpwiedzi skokowej
[a, b, K] = sopdt(y, t, tau, Tp);

%% Wiener
u = [linspace(-45, 45, n); zeros(1, n)];
y = modified_Euler(A, C, alpha_1, alpha_2, F_10, F_D0, h_20, V_10, V_20, 0:Tp:(n-1)*Tp, n, tau, Tp, u);

wiener = Wiener(u, y, a, b, K, n, tau/Tp, mode, s);
wiener = wiener.fuzzyfication();

%% Hammerstein
hammerstein = Hammerstein(45, 135, F_10, F_D0, h_20, alpha_2, n, mode, s);
hammerstein = hammerstein.fuzzyfication();

%% Wymuszenia
u = enforce(kk);
y = modified_Euler(A, C, alpha_1, alpha_2, F_10, F_D0, h_20, V_10, V_20, t, kk, tau, Tp, u);

%% Testowanie
diff_eq(a, b, u, y, kk, tau/Tp, mode)
hammerstein.model(a, b, u, y, K, kk, tau/Tp);
wiener.model(a, b, u, y, K, kk, tau/Tp);