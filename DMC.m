%% Wartość zadana h_2
y_zad = set_value(kk, h_20);

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