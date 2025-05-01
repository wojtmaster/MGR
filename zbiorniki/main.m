%% MGR
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
[s, s_disturbance] = obiekt.linearization(90, 30);
% [a, b] = obiekt.sopdt(obiekt);
% [fis] = obiekt.fuzzyfication();
obiekt.static_charakteristic();

%% Hammerstein
hammerstein = Hammerstein();
hammerstein.linearFuzzy();
hammerstein.nonlinearFuzzy();

% hammerstein.show_fuzzy(hammerstein.linear_fis, 'liniowe');
% hammerstein.show_fuzzy(hammerstein.nonlinear_fis, 'nieliniowe');

%% Wiener
wiener = Wiener(obiekt, a, b);
wiener.linearFuzzy();
wiener.nonlinearFuzzy();

% wiener.show_fuzzy(wiener.linear_fis, 'liniowe');
% wiener.show_fuzzy(wiener.nonlinear_fis, 'nieliniowe');

%% Testing 
for index = 1:1
    % i = 1;
    % U = [repelem((rand(1, obiekt.kk/250) * 90 - 45), 250); 
    %      zeros(1, obiekt.kk)];

    hammerstein.testLinearModel(U, a, b, obiekt, index);
    % hammerstein.testNonlinearModel(U, a, b, obiekt, index);
    
    % wiener.testLinearModel(U, a, b, obiekt, index);
    % wiener.testNonlinearModel(U, a, b, obiekt, index);
end

%% DMC(N, Nu, D, D_disturbance, lambda)
dmc = DMC(100, 2, 150, 150, 1);
% Pierwsza iteracja jest wspólna dla wszystkich algorytmów DMC
dmc.dynamic_matrix(type);
dmc.past_matrix(type);
dmc.matrix_disturbance(s_disturbance);

% Test
y_zad = [repelem((rand(1, obiekt.kk/250) * 40 - 20), 250)];
y_zad(1:100) = 0;
% u_D = [repelem((rand(1, obiekt.kk/200) * 10 - 5), 200)];
u_D = repelem([-5, -3.75, -2.5, -1.25, 0, 1.25, 2.5, 3.75, 5, 5], 200);

%% Wielowątkowy DMC
pool = gcp();
f.analiticH = parfeval(pool, @dmc.dmc_analiticHammerstein, 5, y_zad, u_D, a, b, hammerstein.linear_fis, obiekt.delay, obiekt.kk, 'linear'); 
f.analiticW = parfeval(pool, @dmc.dmc_analiticWiener, 5, y_zad, u_D, a, b, wiener.linear_fis, obiekt.delay, obiekt.kk, 'linear');
f.numericH = parfeval(pool, @dmc.dmc_numericHammerstein, 5, y_zad, u_D, a, b, hammerstein.linear_fis, obiekt.delay, obiekt.kk, 'linear');
f.numericW = parfeval(pool, @dmc.dmc_numericWiener, 5, y_zad, u_D, a, b, wiener.linear_fis, obiekt.delay, obiekt.kk, 'linear');
f.slH = parfeval(pool, @dmc.dmc_slHammerstein, 5, y_zad, u_D, hammerstein.linear_fis, obiekt.delay, obiekt.kk, obiekt.F_10, obiekt.F_D0, @(x, y) obiekt.linearization(x, y), 'linear'); 
f.slW = parfeval(pool, @dmc.dmc_slWiener, 5, y_zad, u_D, wiener.linear_fis, obiekt.delay, obiekt.kk, obiekt.F_10, obiekt.F_D0, @(x, y) obiekt.linearization(x, y), 'linear');
f.nplH = parfeval(pool, @dmc.dmc_nplHammerstein, 5, y_zad, u_D, hammerstein.linear_fis, obiekt.delay, obiekt.kk, obiekt.F_10, obiekt.F_D0, @(x, y) obiekt.linearization(x, y), @(x, y) obiekt.rk4(x, y), 'linear');
f.nplW = parfeval(pool, @dmc.dmc_nplWiener, 5, y_zad, u_D, wiener.linear_fis, obiekt.delay, obiekt.kk, obiekt.F_10, obiekt.F_D0, @(x, y) obiekt.linearization(x, y), @(x, y) obiekt.rk4(x, y), 'linear');
f.fuzzyH = parfeval(pool, @dmc.dmc_fuzzyHammerstein, 5, y_zad, u_D, a, b, hammerstein.linear_fis, obiekt.delay, obiekt.kk, fis, obiekt.F_10, obiekt.F_D0, @(x, y) obiekt.linearization(x, y), 'linear');
f.fuzzyW = parfeval(pool, @dmc.dmc_fuzzyWiener, 5, y_zad, u_D, a, b, wiener.linear_fis, obiekt.delay, obiekt.kk, fis, obiekt.F_10, obiekt.F_D0, @(x, y) obiekt.linearization(x, y), 'linear');

f.analiticHammerstein = parfeval(pool, @dmc.dmc_analiticHammerstein, 5, y_zad, u_D, a, b, hammerstein.nonlinear_fis, obiekt.delay, obiekt.kk, 'nonlinear'); 
f.analiticWiener = parfeval(pool, @dmc.dmc_analiticWiener, 5, y_zad, u_D, a, b, wiener.nonlinear_fis, obiekt.delay, obiekt.kk, 'nonlinear');
f.numericHammerstein = parfeval(pool, @dmc.dmc_numericHammerstein, 5, y_zad, u_D, a, b, hammerstein.nonlinear_fis, obiekt.delay, obiekt.kk, 'nonlinear');
f.numericWiener = parfeval(pool, @dmc.dmc_numericWiener, 5, y_zad, u_D, a, b, wiener.nonlinear_fis, obiekt.delay, obiekt.kk, 'nonlinear');
f.slHammerstein = parfeval(pool, @dmc.dmc_slHammerstein, 5, y_zad, u_D, hammerstein.nonlinear_fis, obiekt.delay, obiekt.kk, obiekt.F_10, obiekt.F_D0, @(x, y) obiekt.linearization(x, y), 'nonlinear'); 
f.slWiener = parfeval(pool, @dmc.dmc_slWiener, 5, y_zad, u_D, wiener.nonlinear_fis, obiekt.delay, obiekt.kk, obiekt.F_10, obiekt.F_D0, @(x, y) obiekt.linearization(x, y), 'nonlinear');
f.nplHammerstein = parfeval(pool, @dmc.dmc_nplHammerstein, 5, y_zad, u_D, hammerstein.nonlinear_fis, obiekt.delay, obiekt.kk, obiekt.F_10, obiekt.F_D0, @(x, y) obiekt.linearization(x, y), @(x, y) obiekt.rk4(x, y), 'nonlinear');
f.nplWiener = parfeval(pool, @dmc.dmc_nplWiener, 5, y_zad, u_D, wiener.nonlinear_fis, obiekt.delay, obiekt.kk, obiekt.F_10, obiekt.F_D0, @(x, y) obiekt.linearization(x, y), @(x, y) obiekt.rk4(x, y), 'nonlinear');
f.fuzzyHammerstein = parfeval(pool, @dmc.dmc_fuzzyHammerstein, 5, y_zad, u_D, a, b, hammerstein.nonlinear_fis, obiekt.delay, obiekt.kk, fis, obiekt.F_10, obiekt.F_D0, @(x, y) obiekt.linearization(x, y), 'nonlinear');
f.fuzzyWiener = parfeval(pool, @dmc.dmc_fuzzyWiener, 5, y_zad, u_D, a, b, wiener.nonlinear_fis, obiekt.delay, obiekt.kk, fis, obiekt.F_10, obiekt.F_D0, @(x, y) obiekt.linearization(x, y), 'nonlinear');

fprintf("\nFunkcje uruchomione asynchronicznie...\n\n");

[y.analiticH, u.analiticH, E.analiticH, E.u_analiticH, E.y_analiticH] = fetchOutputs(f.analiticH);
% dmc.show_result(y.analiticH, y_zad, u.analiticH, E.analiticH, E.u, E.y, obiekt.kk, 'analitic', 'Hammerstein');

[y.analiticW, u.analiticW, E.analiticW, E.u_analiticW, E.y_analiticW] = fetchOutputs(f.analiticW);
% dmc.show_result(y.analiticW, y_zad, u.analiticW, E.analiticW, E.u, E.y, obiekt.kk, 'analitic', 'Wiener');

[y.numericH, u.numericH, E.numericH, E.u_numericH, E.y_numericH] = fetchOutputs(f.numericH);
% dmc.show_result(y.numericH, y_zad, u.numericH, E.numericH, E.u, E.y, obiekt.kk, 'numeric', 'Hammerstein');

[y.numericW, u.numericW, E.numericW, E.u_numericW, E.y_numericW] = fetchOutputs(f.numericW);
% dmc.show_result(y.numericW, y_zad, u.numericW, E.numericW, E.u, E.y, obiekt.kk, 'numeric', 'Wiener');

[y.slH, u.slH, E.slH, E.u_slH, E.y_slH] = fetchOutputs(f.slH);
% dmc.show_result(y.slH, y_zad, u.slH, E.slH, E.u, E.y, obiekt.kk, 'SL', 'Hammerstein');

[y.slW, u.slW, E.slW, E.u_slW, E.y_slW] = fetchOutputs(f.slW);
% dmc.show_result(y.slW, y_zad, u.slW, E.slW, E.u, E.y, obiekt.kk, 'SL', 'Wiener');

[y.nplH, u.nplH, E.nplH, E.u_nplH, E.y_nplH] = fetchOutputs(f.nplH);
% dmc.show_result(y.nplH, y_zad, u.nplH, E.nplH, E.u, E.y, obiekt.kk, 'NPL', 'Hammerstein');

[y.nplW, u.nplW, E.nplW, E.u_nplW, E.y_nplW] = fetchOutputs(f.nplW);
% dmc.show_result(y.nplW, y_zad, u.nplW, E.nplW, E.u, E.y, obiekt.kk, 'NPL', 'Wiener');

[y.fuzzyH, u.fuzzyH, E.fuzzyH, E.u_fuzzyH, E.y_fuzzyH] = fetchOutputs(f.fuzzyH);
% dmc.show_result(y.fuzzyH, y_zad, u.fuzzyH, E.fuzzyH, E.u, E.y, obiekt.kk, 'fuzzy', 'Hammerstein');

[y.fuzzyW, u.fuzzyW, E.fuzzyW, E.u_fuzzyW, E.y_fuzzyW] = fetchOutputs(f.fuzzyW);
% dmc.show_result(y.fuzzyW, y_zad, u.fuzzyW, E.fuzzyW, E.u, E.y, obiekt.kk, 'fuzzy', 'Wiener');

[y.analiticHammerstein, u.analiticHammerstein, E.analiticHammerstein, E.u_analiticHammerstein, E.y_analiticHammerstein] = fetchOutputs(f.analiticHammerstein);
% dmc.show_result(y.analiticHammerstein, y_zad, u.analiticHammerstein, E.analiticHammerstein, E.u, E.y, obiekt.kk, 'analitic', 'Hammerstein');

[y.analiticWiener, u.analiticWiener, E.analiticWiener, E.u_analiticWiener, E.y_analiticWiener] = fetchOutputs(f.analiticWiener);
% dmc.show_result(y.analiticWiener, y_zad, u.analiticWiener, E.analiticWiener, E.u, E.y, obiekt.kk, 'analitic', 'Wiener');

[y.numericHammerstein, u.numericHammerstein, E.numericHammerstein, E.u_numericHammerstein, E.y_numericHammerstein] = fetchOutputs(f.numericHammerstein);
% dmc.show_result(y.numericHammerstein, y_zad, u.numericHammerstein, E.numericHammerstein, E.u, E.y, obiekt.kk, 'numeric', 'Hammerstein');

[y.numericWiener, u.numericWiener, E.numericWiener, E.u_numericWiener, E.y_numericWiener] = fetchOutputs(f.numericWiener);
% dmc.show_result(y.numericWiener, y_zad, u.numericWiener, E.numericWiener, E.u, E.y, obiekt.kk, 'numeric', 'Wiener');

[y.slHammerstein, u.slHammerstein, E.slHammerstein, E.u_slHammerstein, E.y_slHammerstein] = fetchOutputs(f.slHammerstein);
% dmc.show_result(y.slHammerstein, y_zad, u.slHammerstein, E.slHammerstein, E.u, E.y, obiekt.kk, 'SL',' Hammerstein');

[y.slWiener, u.slWiener, E.slWiener, E.u_slWiener, E.y_slWiener] = fetchOutputs(f.slWiener);
% dmc.show_result(y.slWiener, y_zad, u.slWiener, E.slWiener, E.u, E.y, obiekt.kk, 'SL', 'Wiener');

[y.nplHammerstein, u.nplHammerstein, E.nplHammerstein, E.u_nplHammerstein, E.y_nplHammerstein] = fetchOutputs(f.nplHammerstein);
% dmc.show_result(y.nplHammerstein, y_zad, u.nplHammerstein, E.nplHammerstein, E.u, E.y, obiekt.kk, 'NPL', 'Hammerstein');

[y.nplWiener, u.nplWiener, E.nplWiener, E.u_nplWiener, E.y_nplWiener] = fetchOutputs(f.nplWiener);
% dmc.show_result(y.nplWiener, y_zad, u.nplWiener, E.nplWiener, E.u, E.y, obiekt.kk, 'NPL', 'Wiener');

[y.fuzzyHammerstein, u.fuzzyHammerstein, E.fuzzyHammerstein, E.u_fuzzyHammerstein, E.y_fuzzyHammerstein] = fetchOutputs(f.fuzzyHammerstein);
% dmc.show_result(y.fuzzyHammerstein, y_zad, u.fuzzyHammerstein, E.fuzzyHammerstein, E.u, E.y, obiekt.kk, 'fuzzy', 'Hammerstein');

[y.fuzzyWiener, u.fuzzyWiener, E.fuzzyWiener, E.u_fuzzyWiener, E.y_fuzzyWiener] = fetchOutputs(f.fuzzyWiener);
% dmc.show_result(y.fuzzyWiener, y_zad, u.fuzzyWiener, E.fuzzyWiener, E.u, E.y, obiekt.kk, 'fuzzy', 'Wiener');

% delete(pool);

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
stairs(0:length(range)-1, u.(fields{index})(range), 'r-', 'LineWidth', 1.4);
stairs(0:length(range)-1, u.(fields{index+4})(1,(range)), 'g-', 'LineWidth', 1.4);
stairs(0:length(range)-1, u.(fields{index+10})(range), 'b-', 'LineWidth', 1);
stairs(0:length(range)-1, u.(fields{index+14})(1,(range)), 'm-', 'LineWidth', 1);
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
fields = fieldnames(E);
for i = 1:numel(fields)
    disp([fields{i}, ' = ', num2str(E.(fields{i}))]);
end

%%
% Zakresy F1 i FD
F_1 = 0:0.1:45;
F_D = 0:0.1:15;

% Stała alpha_2
alpha_2 = 10; % <- podmień na odpowiednią wartość

% Tworzenie siatki wartości
[F1, FD] = meshgrid(F_1, F_D);

% Obliczanie wartości y
y = ((F1 + FD) / alpha_2).^2;

% Rysowanie wykresu 3D
figure;
surf(F1, FD, y);
xlabel('F_1');
ylabel('F_D');
zlabel('y');
title('Wykres 3D funkcji y = ((F_1 + F_D)/\alpha_2)^2');
shading interp; % (opcjonalnie wygładza kolory)
