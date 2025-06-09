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
% y_zad(1:250) = 0;
% y_zad(251:2000) = 10;
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

[y_analitic, u_analitic, E_analitic, ~, ~] = dmc.dmc_analitic(y_zad, u_D, a, b, obiekt.delay, obiekt.kk);
% [y_no, u_no, E_no, ~, ~] = dmc.dmc_no(y_zad, u_D, a, b, obiekt.delay, obiekt.kk);

[y_sl, u_sl, E_sl, ~, ~] = dmc.dmc_slWiener(y_zad, u_D, a, b, wiener.nonlinear_fis, obiekt.delay, obiekt.kk, 'nonlinear');
[y_npl, u_npl, E_npl, ~, ~] = dmc.dmc_nplWiener(y_zad, u_D, a, b, wiener.nonlinear_fis, obiekt.delay, obiekt.kk, 'nonlinear');
% [y_sl, u_sl, E_sl, ~, ~] = dmc.dmc_slHammerstein(y_zad, u_D, a, b, hammerstein.nonlinear_fis, obiekt.delay, obiekt.kk, 'nonlinear');
% [y_npl, u_npl, E_npl, ~, ~] = dmc.dmc_nplHammerstein(y_zad, u_D, a, b, hammerstein.nonlinear_fis, obiekt.delay, obiekt.kk, 'nonlinear');

% [y_fuzzy, u_fuzzy, E_fuzzy, ~, ~] = dmc.dmc_fuzzy(y_zad, u_D, a, b, obiekt.delay, obiekt.kk, obiekt);
 
t = 0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp;

figure();
plot(t, y_zad, 'k');
hold on;
plot(t, y_analitic, 'b');
plot(t, y_sl, 'r');
plot(t, y_npl, 'g');
% plot(t, y_no, 'c');
% plot(t, y_fuzzy, 'm');
legend('y_zad', 'y_{analitic}', 'y_{sl}', 'y_{npl}', 'y_{no}', 'y_{fuzzy}');
grid on;

% figure();
% stairs(t, u_analitic(1,:), 'b');
% hold on;
% % stairs(t, y_numeric, 'c');
% stairs(t, u_sl(1,:), 'r');
% stairs(t, u_npl(1,:), 'g');
% stairs(t, u_fuzzy(1,:), 'm');
% legend('u_{analitic}', 'u_{sl}', 'u_{npl}', 'u_{fuzzy}');
% grid on;

%% Wielowątkowy DMC
pool = gcp();
f.analiticH = parfeval(pool, @dmc.dmc_analiticHammerstein, 5, y_zad, u_D, a, b, hammerstein.linear_fis, obiekt.delay, obiekt.kk, 'linear'); 
f.analiticW = parfeval(pool, @dmc.dmc_analiticWiener, 5, y_zad, u_D, a, b, wiener.linear_fis, obiekt.delay, obiekt.kk, 'linear');
f.numericH = parfeval(pool, @dmc.dmc_numericHammerstein, 5, y_zad, u_D, a, b, hammerstein.linear_fis, obiekt.delay, obiekt.kk, 'linear');
f.numericW = parfeval(pool, @dmc.dmc_numericWiener, 5, y_zad, u_D, a, b, wiener.linear_fis, obiekt.delay, obiekt.kk, 'linear');
f.slH = parfeval(pool, @dmc.dmc_slHammerstein, 5, y_zad, u_D, hammerstein.linear_fis, obiekt.delay, obiekt.kk, obiekt.F_10, obiekt.F_D0, obiekt, 'linear'); 
f.slW = parfeval(pool, @dmc.dmc_slWiener, 5, y_zad, u_D, wiener.linear_fis, obiekt.delay, obiekt.kk, obiekt.F_10, obiekt.F_D0, obiekt, 'linear');
f.nplH = parfeval(pool, @dmc.dmc_nplHammerstein, 5, y_zad, u_D, hammerstein.linear_fis, obiekt.delay, obiekt.kk, obiekt.F_10, obiekt.F_D0, obiekt, 'linear');
f.nplW = parfeval(pool, @dmc.dmc_nplWiener, 5, y_zad, u_D, wiener.linear_fis, obiekt.delay, obiekt.kk, obiekt.F_10, obiekt.F_D0, obiekt, 'linear');
f.fuzzyH = parfeval(pool, @dmc.dmc_fuzzyHammerstein, 5, y_zad, u_D, a, b, hammerstein.linear_fis, obiekt.delay, obiekt.kk, obiekt.F_10, obiekt.F_D0, obiekt.h_20, obiekt, 'linear');
f.fuzzyW = parfeval(pool, @dmc.dmc_fuzzyWiener, 5, y_zad, u_D, a, b, wiener.linear_fis, obiekt.delay, obiekt.kk, obiekt.F_10, obiekt.F_D0, obiekt.h_20, obiekt, 'linear');

f.analiticHammerstein = parfeval(pool, @dmc.dmc_analiticHammerstein, 5, y_zad, u_D, a, b, hammerstein.nonlinear_fis, obiekt.delay, obiekt.kk, 'nonlinear'); 
f.analiticWiener = parfeval(pool, @dmc.dmc_analiticWiener, 5, y_zad, u_D, a, b, wiener.nonlinear_fis, obiekt.delay, obiekt.kk, 'nonlinear');
f.numericHammerstein = parfeval(pool, @dmc.dmc_numericHammerstein, 5, y_zad, u_D, a, b, hammerstein.nonlinear_fis, obiekt.delay, obiekt.kk, 'nonlinear');
f.numericWiener = parfeval(pool, @dmc.dmc_numericWiener, 5, y_zad, u_D, a, b, wiener.nonlinear_fis, obiekt.delay, obiekt.kk, 'nonlinear');
f.slHammerstein = parfeval(pool, @dmc.dmc_slHammerstein, 5, y_zad, u_D, hammerstein.nonlinear_fis, obiekt.delay, obiekt.kk, obiekt.F_10, obiekt.F_D0, obiekt, 'nonlinear'); 
f.slWiener = parfeval(pool, @dmc.dmc_slWiener, 5, y_zad, u_D, wiener.nonlinear_fis, obiekt.delay, obiekt.kk, obiekt.F_10, obiekt.F_D0, obiekt, 'nonlinear');
f.nplHammerstein = parfeval(pool, @dmc.dmc_nplHammerstein, 5, y_zad, u_D, hammerstein.nonlinear_fis, obiekt.delay, obiekt.kk, obiekt.F_10, obiekt.F_D0, obiekt, 'nonlinear');
f.nplWiener = parfeval(pool, @dmc.dmc_nplWiener, 5, y_zad, u_D, wiener.nonlinear_fis, obiekt.delay, obiekt.kk, obiekt.F_10, obiekt.F_D0, obiekt, 'nonlinear');
f.fuzzyHammerstein = parfeval(pool, @dmc.dmc_fuzzyHammerstein, 5, y_zad, u_D, a, b, hammerstein.nonlinear_fis, obiekt.delay, obiekt.kk, obiekt.F_10, obiekt.F_D0, obiekt.h_20, obiekt, 'nonlinear');
f.fuzzyWiener = parfeval(pool, @dmc.dmc_fuzzyWiener, 5, y_zad, u_D, a, b, wiener.nonlinear_fis, obiekt.delay, obiekt.kk, obiekt.F_10, obiekt.F_D0, obiekt.h_20, obiekt, 'nonlinear');

fprintf("\nFunkcje uruchomione asynchronicznie...\n\n");

[y.analiticH, u_sl.analiticH, E_sl.analiticH, E_sl.u_analiticH, E_sl.y_analiticH] = fetchOutputs(f.analiticH);
% dmc.show_result(y.analiticH, y_zad, u.analiticH, E.analiticH, E.u_analiticH, E.y_analiticH, obiekt.kk, 'analitic', 'Hammerstein');

[y.analiticW, u_sl.analiticW, E_sl.analiticW, E_sl.u_analiticW, E_sl.y_analiticW] = fetchOutputs(f.analiticW);
% dmc.show_result(y.analiticW, y_zad, u.analiticW, E.analiticW, E.u_analiticW, E.y_analiticW, obiekt.kk, 'analitic', 'Wiener');

[y.numericH, u_sl.numericH, E_sl.numericH, E_sl.u_numericH, E_sl.y_numericH] = fetchOutputs(f.numericH);
% dmc.show_result(y.numericH, y_zad, u.numericH, E.numericH, E.u_numericH, E.y_numericH, obiekt.kk, 'numeric', 'Hammerstein');

[y.numericW, u_sl.numericW, E_sl.numericW, E_sl.u_numericW, E_sl.y_numericW] = fetchOutputs(f.numericW);
% dmc.show_result(y.numericW, y_zad, u.numericW, E.numericW, E.u_numericW, E.y_numericW, obiekt.kk, 'numeric', 'Wiener');

[y.slH, u_sl.slH, E_sl.slH, E_sl.u_slH, E_sl.y_slH] = fetchOutputs(f.slH);
% dmc.show_result(y.slH, y_zad, u.slH, E.slH, E.u_slH, E.y_slH, obiekt.kk, 'SL', 'Hammerstein');

[y.slW, u_sl.slW, E_sl.slW, E_sl.u_slW, E_sl.y_slW] = fetchOutputs(f.slW);
% dmc.show_result(y.slW, y_zad, u.slW, E.slW, E.u_slW, E.y_slW, obiekt.kk, 'SL', 'Wiener');

[y.nplH, u_sl.nplH, E_sl.nplH, E_sl.u_nplH, E_sl.y_nplH] = fetchOutputs(f.nplH);
% dmc.show_result(y.nplH, y_zad, u.nplH, E.nplH, E.u_nplH, E.y_nplH, obiekt.kk, 'NPL', 'Hammerstein');

[y.nplW, u_sl.nplW, E_sl.nplW, E_sl.u_nplW, E_sl.y_nplW] = fetchOutputs(f.nplW);
% dmc.show_result(y.nplW, y_zad, u.nplW, E.nplW, E.u_nplW, E.y_nplW, obiekt.kk, 'NPL', 'Wiener');

[y.fuzzyH, u_sl.fuzzyH, E_sl.fuzzyH, E_sl.u_fuzzyH, E_sl.y_fuzzyH] = fetchOutputs(f.fuzzyH);
% dmc.show_result(y.fuzzyH, y_zad, u.fuzzyH, E.fuzzyH, E.u_fuzzyH, E.y_fuzzyH, obiekt.kk, 'fuzzy', 'Hammerstein');

[y.fuzzyW, u_sl.fuzzyW, E_sl.fuzzyW, E_sl.u_fuzzyW, E_sl.y_fuzzyW] = fetchOutputs(f.fuzzyW);
% dmc.show_result(y.fuzzyW, y_zad, u.fuzzyW, E.fuzzyW, E.u_fuzzyW, E.y_fuzzyW, obiekt.kk, 'fuzzy', 'Wiener');

[y.analiticHammerstein, u_sl.analiticHammerstein, E_sl.analiticHammerstein, E_sl.u_analiticHammerstein, E_sl.y_analiticHammerstein] = fetchOutputs(f.analiticHammerstein);
% dmc.show_result(y.analiticHammerstein, y_zad, u.analiticHammerstein, E.analiticHammerstein, E.u_analiticHammerstein, E.y_analiticHammerstein, obiekt.kk, 'analitic', 'Hammerstein');

[y.analiticWiener, u_sl.analiticWiener, E_sl.analiticWiener, E_sl.u_analiticWiener, E_sl.y_analiticWiener] = fetchOutputs(f.analiticWiener);
% dmc.show_result(y.analiticWiener, y_zad, u.analiticWiener, E.analiticWiener, E.u_analiticWiener, E.y_analiticWiener, obiekt.kk, 'analitic', 'Wiener');

[y.numericHammerstein, u_sl.numericHammerstein, E_sl.numericHammerstein, E_sl.u_numericHammerstein, E_sl.y_numericHammerstein] = fetchOutputs(f.numericHammerstein);
% dmc.show_result(y.numericHammerstein, y_zad, u.numericHammerstein, E.numericHammerstein, E.u_numericHammerstein, E.y_numericHammerstein, obiekt.kk, 'numeric', 'Hammerstein');

[y.numericWiener, u_sl.numericWiener, E_sl.numericWiener, E_sl.u_numericWiener, E_sl.y_numericWiener] = fetchOutputs(f.numericWiener);
% dmc.show_result(y.numericWiener, y_zad, u.numericWiener, E.numericWiener, E.u_numericWiener, E.y_numericWiener, obiekt.kk, 'numeric', 'Wiener');

[y.slHammerstein, u_sl.slHammerstein, E_sl.slHammerstein, E_sl.u_slHammerstein, E_sl.y_slHammerstein] = fetchOutputs(f.slHammerstein);
% dmc.show_result(y.slHammerstein, y_zad, u.slHammerstein, E.slHammerstein, E.u_slHammerstein, E.y_slHammerstein, obiekt.kk, 'SL',' Hammerstein');

[y.slWiener, u_sl.slWiener, E_sl.slWiener, E_sl.u_slWiener, E_sl.y_slWiener] = fetchOutputs(f.slWiener);
% % dmc.show_result(y.slWiener, y_zad, u.slWiener, E.slWiener, E.u_slWiener, E.y_slWiener, obiekt.kk, 'SL', 'Wiener');

[y.nplHammerstein, u_sl.nplHammerstein, E_sl.nplHammerstein, E_sl.u_nplHammerstein, E_sl.y_nplHammerstein] = fetchOutputs(f.nplHammerstein);
% dmc.show_result(y.nplHammerstein, y_zad, u.nplHammerstein, E.nplHammerstein, E.u_nplHammerstein, E.y_nplHammerstein, obiekt.kk, 'NPL', 'Hammerstein');

[y.nplWiener, u_sl.nplWiener, E_sl.nplWiener, E_sl.u_nplWiener, E_sl.y_nplWiener] = fetchOutputs(f.nplWiener);
% dmc.show_result(y.nplWiener, y_zad, u.nplWiener, E.nplWiener, E.u_nplWiener, E.y_nplWiener, obiekt.kk, 'NPL', 'Wiener');

[y.fuzzyHammerstein, u_sl.fuzzyHammerstein, E_sl.fuzzyHammerstein, E_sl.u_fuzzyHammerstein, E_sl.y_fuzzyHammerstein] = fetchOutputs(f.fuzzyHammerstein);
% dmc.show_result(y.fuzzyHammerstein, y_zad, u.fuzzyHammerstein, E.fuzzyHammerstein, E.u_fuzzyHammerstein, E.y_fuzzyHammerstein, obiekt.kk, 'fuzzy', 'Hammerstein');

[y.fuzzyWiener, u_sl.fuzzyWiener, E_sl.fuzzyWiener, E_sl.u_fuzzyWiener, E_sl.y_fuzzyWiener] = fetchOutputs(f.fuzzyWiener);
% dmc.show_result(y.fuzzyWiener, y_zad, u.fuzzyWiener, E.fuzzyWiener, E.u_fuzzyWiener, E.y_fuzzyWiener, obiekt.kk, 'fuzzy', 'Wiener');

% delete(pool);

%%
t = 0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp;
figure;
plot(t, y_zad, 'k-');
hold on;
plot(t, y.analiticHammerstein, 'b-');
plot(t, y.numericHammerstein, 'r-');
plot(t, y.slHammerstein, 'g-');
plot(t, y.nplHammerstein, 'm-');
plot(t, y.fuzzyHammerstein, 'c-');
legend('y_{zad}', 'y_{analiticHammerstein}', ...
    'y_{numericHammerstein}', 'y_{slHammerstein}', ...
    'y_{nplHammerstein}', 'y_{fuzzyHammerstein}', 'Location', 'best');
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
plot(t, y.fuzzyWiener, 'c-');
legend('y_{zad}', 'y_{analiticWiener}', ...
    'y_{numericWiener}', 'y_{slWiener}', ...
    'y_{nplWiener}', 'y_{fuzzyWiener}', 'Location', 'best');
title(sprintf('Sygnał wyjściowy y(k) - Wiener'));
xlabel('k');
ylabel('y(k)');
grid on;

%% Prezentacja wyników
fields = fieldnames(y);
%%
range = 1:obiekt.kk;
type = 'numeric';
index = 1;

figure;
hold on;
stairs(0:length(range)-1, y.(fields{index})(range), 'r-', 'LineWidth', 1.4);
stairs(0:length(range)-1, y.(fields{index+4})(range), 'g-', 'LineWidth', 1.4);
stairs(0:length(range)-1, y.(fields{index+10})(range), 'b-', 'LineWidth', 1);
stairs(0:length(range)-1, y.(fields{index+14})(range), 'm-', 'LineWidth', 1);
stairs(0:length(range)-1, y_zad(range), 'k-', 'LineWidth', 1);
hold off;
legend('Hammerstein (następniki liniowe)', ...
    'Wiener (następniki liniowe)', ...
    'Hammerstein (następniki nieliniowe)', ...
    'Wiener (następniki nieliniowe)', ...,
    'y_{zad}', 'Location', 'best');
title(sprintf('Sygnał wyjściowy y(k) \t DMC - %s \t %s', type));
xlabel('k');
ylabel('y(k)');
grid on;
% saveas(gcf, sprintf('D:/EiTI/MGR/raporty/raport_MGR/pictures/DMC_%s_y.png', s));  % Zapisuje jako plik PNG

figure;
hold on;
stairs(0:length(range)-1, u_sl.(fields{index})(range), 'r-', 'LineWidth', 1.4);
stairs(0:length(range)-1, u_sl.(fields{index+4})(1,(range)), 'g-', 'LineWidth', 1.4);
stairs(0:length(range)-1, u_sl.(fields{index+10})(range), 'b-', 'LineWidth', 1);
stairs(0:length(range)-1, u_sl.(fields{index+14})(1,(range)), 'm-', 'LineWidth', 1);
hold off;
legend('Hammerstein (następniki liniowe)', ...
    'Wiener (następniki liniowe)', ...
    'Hammerstein (następniki nieliniowe)', ...
    'Wiener (następniki nieliniowe)', ...,
    'Location', 'best');
title(sprintf('Sygnał sterujący u(k) \t DMC - %s \t %s', type));
xlabel('k');
ylabel('u(k)');
grid on;
% saveas(gcf, sprintf('D:/EiTI/MGR/raporty/raport_MGR/pictures/DMC_%s_u.png', s));  % Zapisuje jako plik PNG

%% Errors
fields = fieldnames(E_sl);
for i = 1:numel(fields)
    disp([fields{i}, ' = ', num2str(E_sl.(fields{i}))]);
end