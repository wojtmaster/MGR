% %% MGR
% Autor: Wojciech Rogalski
% Data: 15.12.2024r.
% Tytuł: Porównanie modeli Hammersteina i Winera w regulacji rozmytej
clear;
close;
clc;
set(groot, 'DefaultTextFontName', 'Arial');
set(groot, 'DefaultAxesFontName', 'Arial');
set(groot, 'DefaultTextFontSize', 12);
set(groot, 'DefaultAxesFontSize', 10);

%% Inicjalizacja obiektu oraz linearyzacja
obiekt = Obiekt();
obiekt.linearization(90, 30);
[a, b, ~, ~] = obiekt.sopdt();
% obiekt.static_charakteristic();

%% Hammerstein
hammerstein = Hammerstein(obiekt);
hammerstein.linearFuzzy();
hammerstein.nonlinearFuzzy();

% hammerstein.show_fuzzy(hammerstein.linear_fis, 'liniowe');
% hammerstein.show_fuzzy(hammerstein.nonlinear_fis, 'nieliniowe');

%% Wiener
wiener = Wiener(a, b, obiekt);
wiener.linearFuzzy();
wiener.nonlinearFuzzy();

% wiener.show_fuzzy(wiener.linear_fis, 'liniowe');
% wiener.show_fuzzy(wiener.nonlinear_fis, 'nieliniowe');

%% Testing 
for index = 1:1
    % i = 1;
    U = [repelem((rand(1, obiekt.kk/200) * 80 - 35), 200); 
         repelem((rand(1, obiekt.kk/250) * 30 - 15), 250)];

    hammerstein.testLinearModel(U, a, b, obiekt, index);
    hammerstein.testNonlinearModel(U, a, b, obiekt, index);

    wiener.testLinearModel(U, a, b, obiekt, index);
    wiener.testNonlinearModel(U, a, b, obiekt, index);
end

%% DMC(N, Nu, D, D_disturbance, lambda)
close all;
y_zad = [repelem((rand(1, obiekt.kk/400) * 20 - 10), 400)];
y_zad(1:100) = 0;
u_D = [repelem((rand(1, obiekt.kk/250) * 10 - 5), 250)];
u_D(1:obiekt.kk) = 0;

obiekt.linearization(90, 30);
[a, b, s, s_D] = obiekt.sopdt();
s(1) = [];
s_D(1) = [];
s(obiekt.dynamic_horizont) = s(end);
s_D(obiekt.dynamic_horizont) = s_D(end);

dmc = DMC(50, 25, 150, 150, 50);
dmc.dynamic_matrix(s); 
dmc.past_matrix(s);
dmc.matrix_disturbance(s_D);

% [y_analitic, u_analitic, E_analitic, ~, ~] = dmc.dmc_analiticWiener(y_zad, u_D, a, b, obiekt);
% [y_no, u_no, E_no, ~, ~] = dmc.dmc_noWiener(y_zad, u_D, a, b, obiekt);
% [y_sl, u_sl, E_sl, ~, ~] = dmc.dmc_slWiener(y_zad, u_D, a, b, wiener.nonlinear_fis, obiekt, 'nonlinear');
% [y_npl, u_npl, E_npl, ~, ~] = dmc.dmc_nplWiener(y_zad, u_D, a, b, wiener.nonlinear_fis, obiekt, 'nonlinear');
% [y_fuzzy, u_fuzzy, E_fuzzy, ~, ~] = dmc.dmc_fuzzyWiener(y_zad, u_D, a, b, obiekt);

% [y_analitic, u_analitic, E_analitic, ~, ~] = dmc.dmc_analiticHammerstein(y_zad, u_D, a, b, obiekt);
% [y_no, u_no, E_no, ~, ~] = dmc.dmc_noHammerstein(y_zad, u_D, a, b, obiekt);
% [y_sl, u_sl, E_sl, ~, ~] = dmc.dmc_slHammerstein(y_zad, u_D, a, b, hammerstein.nonlinear_fis, obiekt, 'nonlinear');
% [y_npl, u_npl, E_npl, ~, ~] = dmc.dmc_nplHammerstein(y_zad, u_D, a, b, hammerstein.nonlinear_fis, obiekt, 'nonlinear');
% [y_fuzzy, u_fuzzy, E_fuzzy, ~, ~] = dmc.dmc_fuzzyHammerstein(y_zad, u_D, a, b, obiekt);
% 
% t = 0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp;
% 
% figure();
% plot(t, y_zad, 'k');
% hold on;
% plot(t, y_analitic, 'b');
% plot(t, y_sl, 'r');
% plot(t, y_npl, 'g');
% plot(t, y_no, 'c');
% plot(t, y_fuzzy, 'm');
% legend('y_zad', 'y_{analitic}', 'y_{sl}', 'y_{npl}', 'y_{no}', 'y_{fuzzy}');
% grid on;
% 
% figure();
% stairs(t, u_analitic(1,:), 'b');
% hold on;
% % stairs(t, u_numeric(1,:), 'c');
% stairs(t, u_sl(1,:), 'r');
% stairs(t, u_npl(1,:), 'g');
% stairs(t, u_fuzzy(1,:), 'm');
% legend('u_{analitic}', 'u_{sl}', 'u_{npl}', 'u_{fuzzy}');
% grid on;

%% Wielowątkowy DMC
pool = gcp();

b(2) = 0;

f.analiticHammerstein = parfeval(pool, @dmc.dmc_analiticHammerstein, 5, y_zad, u_D, a, b, obiekt); 
f.numericHammerstein = parfeval(pool, @dmc.dmc_numericHammerstein, 5, y_zad, u_D, a, b, obiekt);
f.noHammerstein = parfeval(pool, @dmc.dmc_noHammerstein, 5, y_zad, u_D, a, b, obiekt);
f.fuzzyHammerstein = parfeval(pool, @dmc.dmc_fuzzyHammerstein, 5, y_zad, u_D, a, b, obiekt);
f.slH = parfeval(pool, @dmc.dmc_slHammerstein, 5, y_zad, u_D, a, b, hammerstein.linear_fis, obiekt, 'linear'); 
f.slHammerstein = parfeval(pool, @dmc.dmc_slHammerstein, 5, y_zad, u_D, a, b, hammerstein.nonlinear_fis, obiekt, 'nonlinear');
f.nplH = parfeval(pool, @dmc.dmc_nplHammerstein, 5, y_zad, u_D, a, b, hammerstein.linear_fis, obiekt, 'linear');
f.nplHammerstein = parfeval(pool, @dmc.dmc_nplHammerstein, 5, y_zad, u_D, a, b, hammerstein.nonlinear_fis, obiekt, 'nonlinear');

f.analiticWiener = parfeval(pool, @dmc.dmc_analiticWiener, 5, y_zad, u_D, a, b, obiekt); 
f.numericWiener = parfeval(pool, @dmc.dmc_numericWiener, 5, y_zad, u_D, a, b, obiekt);
f.noWiener = parfeval(pool, @dmc.dmc_noWiener, 5, y_zad, u_D, a, b, obiekt);
f.fuzzyWiener = parfeval(pool, @dmc.dmc_fuzzyWiener, 5, y_zad, u_D, a, b, obiekt);
f.slW = parfeval(pool, @dmc.dmc_slWiener, 5, y_zad, u_D, a, b, wiener.linear_fis, obiekt, 'linear');
f.slWiener = parfeval(pool, @dmc.dmc_slWiener, 5, y_zad, u_D, a, b, wiener.nonlinear_fis, obiekt, 'nonlinear');
f.nplW = parfeval(pool, @dmc.dmc_nplWiener, 5, y_zad, u_D, a, b, wiener.linear_fis, obiekt, 'linear');
f.nplWiener = parfeval(pool, @dmc.dmc_nplWiener, 5, y_zad, u_D, a, b, wiener.nonlinear_fis, obiekt, 'nonlinear');

fprintf("\nFunkcje uruchomione asynchronicznie...\n\n");

[y.analiticHammerstein, u.analiticHammerstein, E.analiticHammerstein, E.u_analiticHammerstein, E.y_analiticHammerstein] = fetchOutputs(f.analiticHammerstein);
% dmc.show_result(y.analiticHammerstein, y_zad, u.analiticHammerstein, E.analiticHammerstein, E.u_analiticHammerstein, E.y_analiticHammerstein, obiekt.kk, 'analitic', 'Hammerstein');

[y.analiticWiener, u.analiticWiener, E.analiticWiener, E.u_analiticWiener, E.y_analiticWiener] = fetchOutputs(f.analiticWiener);
% dmc.show_result(y.analiticWiener, y_zad, u.analiticWiener, E.analiticWiener, E.u_analiticWiener, E.y_analiticWiener, obiekt.kk, 'analitic', 'Wiener');

[y.numericHammerstein, u.numericHammerstein, E.numericHammerstein, E.u_numericHammerstein, E.y_numericHammerstein] = fetchOutputs(f.numericHammerstein);
% dmc.show_result(y.numericHammerstein, y_zad, u.numericHammerstein, E.numericHammerstein, E.u_numericHammerstein, E.y_numericHammerstein, obiekt.kk, 'numeric', 'Hammerstein');

[y.numericWiener, u.numericWiener, E.numericWiener, E.u_numericWiener, E.y_numericWiener] = fetchOutputs(f.numericWiener);
% dmc.show_result(y.numericWiener, y_zad, u.numericWiener, E.numericWiener, E.u_numericWiener, E.y_numericWiener, obiekt.kk, 'numeric', 'Wiener');

[y.noHammerstein, u.noHammerstein, E.noHammerstein, E.u_noHammerstein, E.y_noHammerstein] = fetchOutputs(f.noHammerstein);
% dmc.show_result(y.noHammerstein, y_zad, u.noHammerstein, E.noHammerstein, E.u_noHammerstein, E.y_noHammerstein, obiekt.kk, 'no', 'Hammerstein');

[y.noWiener, u.noWiener, E.noWiener, E.u_noWiener, E.y_noWiener] = fetchOutputs(f.noWiener);
% dmc.show_result(y.noWiener, y_zad, u.noWiener, E.noWiener, E.u_noWiener, E.y_noWiener, obiekt.kk, 'no', 'Wiener');

[y.fuzzyHammerstein, u.fuzzyHammerstein, E.fuzzyHammerstein, E.u_fuzzyHammerstein, E.y_fuzzyHammerstein] = fetchOutputs(f.fuzzyHammerstein);
% dmc.show_result(y.fuzzyHammerstein, y_zad, u.fuzzyHammerstein, E.fuzzyHammerstein, E.u_fuzzyHammerstein, E.y_fuzzyHammerstein, obiekt.kk, 'fuzzy', 'Hammerstein');

[y.fuzzyWiener, u.fuzzyWiener, E.fuzzyWiener, E.u_fuzzyWiener, E.y_fuzzyWiener] = fetchOutputs(f.fuzzyWiener);
% dmc.show_result(y.fuzzyWiener, y_zad, u.fuzzyWiener, E.fuzzyWiener, E.u_fuzzyWiener, E.y_fuzzyWiener, obiekt.kk, 'fuzzy', 'Wiener');

[y.slHammerstein, u.slHammerstein, E.slHammerstein, E.u_slHammerstein, E.y_slHammerstein] = fetchOutputs(f.slHammerstein);
% dmc.show_result(y.slHammerstein, y_zad, u.slHammerstein, E.slHammerstein, E.u_slHammerstein, E.y_slHammerstein, obiekt.kk, 'SL',' Hammerstein');

[y.slWiener, u.slWiener, E.slWiener, E.u_slWiener, E.y_slWiener] = fetchOutputs(f.slWiener);
% % dmc.show_result(y.slWiener, y_zad, u.slWiener, E.slWiener, E.u_slWiener, E.y_slWiener, obiekt.kk, 'SL', 'Wiener');

[y.nplHammerstein, u.nplHammerstein, E.nplHammerstein, E.u_nplHammerstein, E.y_nplHammerstein] = fetchOutputs(f.nplHammerstein);
% dmc.show_result(y.nplHammerstein, y_zad, u.nplHammerstein, E.nplHammerstein, E.u_nplHammerstein, E.y_nplHammerstein, obiekt.kk, 'NPL', 'Hammerstein');

[y.nplWiener, u.nplWiener, E.nplWiener, E.u_nplWiener, E.y_nplWiener] = fetchOutputs(f.nplWiener);
% dmc.show_result(y.nplWiener, y_zad, u.nplWiener, E.nplWiener, E.u_nplWiener, E.y_nplWiener, obiekt.kk, 'NPL', 'Wiener');

[y.slH, u.slH, E.slH, E.u_slH, E.y_slH] = fetchOutputs(f.slH);
% dmc.show_result(y.slH, y_zad, u.slH, E.slH, E.u_slH, E.y_slH, obiekt.kk, 'SL', 'Hammerstein');

[y.slW, u.slW, E.slW, E.u_slW, E.y_slW] = fetchOutputs(f.slW);
% dmc.show_result(y.slW, y_zad, u.slW, E.slW, E.u_slW, E.y_slW, obiekt.kk, 'SL', 'Wiener');

[y.nplH, u.nplH, E.nplH, E.u_nplH, E.y_nplH] = fetchOutputs(f.nplH);
% dmc.show_result(y.nplH, y_zad, u.nplH, E.nplH, E.u_nplH, E.y_nplH, obiekt.kk, 'NPL', 'Hammerstein');

[y.nplW, u.nplW, E.nplW, E.u_nplW, E.y_nplW] = fetchOutputs(f.nplW);
% dmc.show_result(y.nplW, y_zad, u.nplW, E.nplW, E.u_nplW, E.y_nplW, obiekt.kk, 'NPL', 'Wiener');

% delete(pool);

%% Prezentacja wyników
t = 0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp;
figure;
plot(t, y_zad, 'k-');
hold on;
plot(t, y.analiticHammerstein, 'b-');
plot(t, y.numericHammerstein, 'r-');
plot(t, y.slHammerstein, 'g-');
plot(t, y.nplHammerstein, 'm-');
plot(t, y.noHammerstein, 'y-');
plot(t, y.fuzzyHammerstein, 'c-');
legend('y_{zad}', 'y_{analitic}', ...
    'y_{numeric}', 'y_{slHammerstein}', ...
    'y_{nplHammerstein}', 'y_{no}', 'y_{fuzzy}', 'Location', 'best');
title(sprintf('Sygnał wyjściowy y(k) - Hammerstein'));
xlabel('k');
ylabel('y(k)');
grid on;

figure;
plot(t, y_zad, 'k-');
hold on;
plot(t, y.analiticWiener, 'b-');
plot(t, y.numericWiener, 'r-');
plot(t, y.slWiener, 'g-');
plot(t, y.nplWiener, 'm-');
plot(t, y.noWiener, 'y-');
plot(t, y.fuzzyWiener, 'c-');
legend('y_{zad}', 'y_{analitic}', ...
    'y_{numeric}', 'y_{slWiener}', ...
    'y_{nplWiener}', 'y_{no}', 'y_{fuzzy}', 'Location', 'best');
title(sprintf('Sygnał wyjściowy y(k) - Wiener'));
xlabel('k');
ylabel('y(k)');
grid on;

figure;
plot(t, y_zad, 'k-');
hold on;
plot(t, y.analiticHammerstein, 'b-');
plot(t, y.numericHammerstein, 'r-');
plot(t, y.slH, 'g-');
plot(t, y.nplH, 'm-');
plot(t, y.noHammerstein, 'y-');
plot(t, y.fuzzyHammerstein, 'c-');
legend('y_{zad}', 'y_{analitic}', ...
    'y_{numeric}', 'y_{slHammerstein}', ...
    'y_{nplHammerstein}', 'y_{no}', 'y_{fuzzy}', 'Location', 'best');
title(sprintf('Sygnał wyjściowy y(k) - Hammerstein'));
xlabel('k');
ylabel('y(k)');
grid on;

figure;
plot(t, y_zad, 'k-');
hold on;
plot(t, y.analiticWiener, 'b-');
plot(t, y.numericWiener, 'r-');
plot(t, y.slW, 'g-');
plot(t, y.nplW, 'm-');
plot(t, y.noWiener, 'y-');
plot(t, y.fuzzyWiener, 'c-');
legend('y_{zad}', 'y_{analitic}', ...
    'y_{numeric}', 'y_{slWiener}', ...
    'y_{nplWiener}', 'y_{no}', 'y_{fuzzy}', 'Location', 'best');
title(sprintf('Sygnał wyjściowy y(k) - Wiener'));
xlabel('k');
ylabel('y(k)');
grid on;