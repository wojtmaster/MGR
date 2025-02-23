%% MGR
% Autor: Wojciech Rogalski
% Data: 15.12.2024r.
% Tytuł: Porównanie modeli Hammersteina i Winera w regulacji rozmytej
clear all;
clc;
set(groot, 'DefaultTextFontName', 'Arial');
set(groot, 'DefaultAxesFontName', 'Arial');
set(groot, 'DefaultTextFontSize', 12);
set(groot, 'DefaultAxesFontSize', 10);

%% Inicjalizacja obiektu oraz linearyzacja
obiekt = Obiekt();
[s, s_disturbance, a, b] = obiekt.linearization(90, 30);
% [fis] = obiekt.fuzzyfication();

%% DMC(N, Nu, D, D_disturbance, lambda)
dmc = DMC(100, 2, 150, 150, 1);
% Pierwsza iteracja jest wspólna dla wszystkich algorytmów DMC
dmc.dynamic_matrix(s);
dmc.past_matrix(s);
dmc.matrix_disturbance(s_disturbance);

% Test
y_zad = [repelem((rand(1, obiekt.kk/250) * 40 - 20), 250)];
y_zad(1:100) = 0;
u_D = [repelem((rand(1, obiekt.kk/200) * 10 - 5), 200)];

pool = gcp();
f.analitic = parfeval(pool, @dmc.dmc_analitic, 3, y_zad, u_D, a, b, obiekt.delay, obiekt.kk);
f.numeric = parfeval(pool, @dmc.dmc_numeric, 3, y_zad, u_D, a, b, obiekt.delay, obiekt.kk);
f.sl = parfeval(pool, @dmc.dmc_sl, 3, y_zad, u_D, obiekt.delay, obiekt.kk, obiekt.F_10, obiekt.F_D0, @(x, y) obiekt.linearization(x, y));
f.SL = parfeval(pool, @dmc.dmc_SL, 3, y_zad, u_D, obiekt.delay, obiekt.kk, obiekt.F_10, obiekt.F_D0, @(x, y) obiekt.linearization(x, y));
f.npl = parfeval(pool, @dmc.dmc_npl, 3, y_zad, u_D, obiekt.delay, obiekt.kk, obiekt.F_10, obiekt.F_D0, @(x, y) obiekt.linearization(x, y), @(x, y) obiekt.rk4(x, y));
f.NPL = parfeval(pool, @dmc.dmc_NPL, 3, y_zad, u_D, obiekt.delay, obiekt.kk, obiekt.F_10, obiekt.F_D0, @(x, y) obiekt.linearization(x, y), @(x, y) obiekt.rk4(x, y));
f.fuzzy = parfeval(pool, @dmc.dmc_fuzzy, 3, y_zad, u_D, a, b, obiekt.delay, obiekt.kk, fis, obiekt.F_10, obiekt.F_D0, @(x, y) obiekt.linearization(x, y));
f.FUZZY = parfeval(pool, @dmc.dmc_FUZZY, 3, y_zad, u_D, a, b, obiekt.delay, obiekt.kk, fis, obiekt.F_10, obiekt.F_D0, @(x, y) obiekt.linearization(x, y));
fprintf("\nFunkcje uruchomione asynchronicznie...\n\n");

[y.analitic, u.analitic, E.analitic] = fetchOutputs(f.analitic);
dmc.show_result(y.analitic, y_zad, u.analitic, E.analitic, obiekt.kk, 'analitic');

[y.numeric, u.numeric, E.numeric] = fetchOutputs(f.numeric);
dmc.show_result(y.numeric, y_zad, u.numeric, E.numeric, obiekt.kk, 'numeric');

[y.sl, u.sl, E.sl] = fetchOutputs(f.sl);
dmc.show_result(y.sl, y_zad, u.sl, E.sl, obiekt.kk, 'sl');

[y.SL, u.SL, E.SL] = fetchOutputs(f.SL);
dmc.show_result(y.SL, y_zad, u.SL, E.SL, obiekt.kk, 'SL');

[y.npl, u.npl, E.npl] = fetchOutputs(f.npl);
dmc.show_result(y.npl, y_zad, u.npl, E.npl, obiekt.kk, 'npl');

[y.NPL, u.NPL, E.NPL] = fetchOutputs(f.NPL);
dmc.show_result(y.NPL, y_zad, u.NPL, E.NPL, obiekt.kk, 'NPL');

[y.fuzzy, u.fuzzy, E.fuzzy] = fetchOutputs(f.fuzzy);
dmc.show_result(y.fuzzy, y_zad, u.fuzzy, E.fuzzy, obiekt.kk, 'fuzzy');

[y.FUZZY, u.FUZZY, E.FUZZY] = fetchOutputs(f.FUZZY);
dmc.show_result(y.FUZZY, y_zad, u.FUZZY, E.FUZZY, obiekt.kk, 'FUZZY');

delete(pool);

%% Hammerstein
hammerstein = Hammerstein();
hammerstein.linearFuzzy();
hammerstein.nonlinearFuzzy();

% hammerstein.show_fuzzy(hammerstein.linear_fis);
% hammerstein.show_fuzzy(hammerstein.nonlinear_fis);

%% Wiener
wiener = Wiener();
wiener.linearFuzzy();
wiener.nonlinearFuzzy();

% wiener.show_fuzzy(wiener.linear_fis);
% wiener.show_fuzzy(wiener.nonlinear_fis);

%% Testing 
U = [repelem((rand(1, obiekt.kk/400) * 90 - 45), 400); zeros(1, obiekt.kk)];

hammerstein.testLinearModel(U, a, b, obiekt.delay, obiekt.kk, obiekt.Tp, @(x, y) obiekt.rk4(x, y));
hammerstein.testNonlinearModel(U, a, b, obiekt.delay, obiekt.kk, obiekt.Tp, @(x, y) obiekt.rk4(x, y));

wiener.testLinearModel(U, a, b, obiekt.delay, obiekt.kk, obiekt.Tp, @(x, y) obiekt.rk4(x, y));
wiener.testNonlinearModel(U, a, b, obiekt.delay, obiekt.kk, obiekt.Tp, @(x, y) obiekt.rk4(x, y));