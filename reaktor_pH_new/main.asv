%% MGR
% Autor: Wojciech Rogalski
% Data: 15.12.2024r.
% Tytuł: Porównanie modeli Hammersteina i Winera w regulacji predykcyjnej
clear all;
clc;
set(groot, 'DefaultTextFontName', 'Arial');
set(groot, 'DefaultAxesFontName', 'Arial');
set(groot, 'DefaultTextFontSize', 12);
set(groot, 'DefaultAxesFontSize', 10);

%% Inicjalizacja obiektu oraz linearyzacja
obiekt = Obiekt();
%  pH_0, h_0, Q_10, Q_20, Q_30
obiekt.linearization(16.6, 0.55, 15.6);
% obiekt.show_staticCharacteristic();

[a_h, b_h, s_h] = obiekt.fopdtModel('h', 1);
[a_pH, b_pH, s_pH] = obiekt.tfestModel();

[a_Wa4, b_Wa4, ~] = obiekt.fopdtModel('W_{a4}', 3);
[a_Wb4, b_Wb4, ~] = obiekt.fopdtModel('W_{b4}', 4);

%% Hammerstein 
hammerstein = Hammerstein(obiekt);
hammerstein.linearFuzzy();
hammerstein.nonlinearFuzzy();
% hammerstein.checkStatic('linear');
% hammerstein.checkStatic('nonlinear');

%% Wiener
wiener = Wiener(a_h, a_Wa4, a_Wb4, b_h, b_Wa4, b_Wb4, obiekt);
wiener.linearFuzzy();
wiener.nonlinearFuzzy();
% wiener.checkStatic('linear');
% wiener.checkStatic('nonlinear');

%% Testy Hammerstein
U = [repelem((rand(1, obiekt.kk/400) * 30 - 15), 400);
    zeros(1, obiekt.kk);
    repelem((rand(1, obiekt.kk/600) * 30 - 15), 600)];

hammerstein.testLinearModel(U, a_h, a_pH, b_h, b_pH, obiekt, 'h', 1);
hammerstein.testNonlinearModel(U, a_h, a_pH, b_h, b_pH, obiekt, 'h', 1);

hammerstein.testLinearModel(U, a_h, a_pH, b_h, b_pH, obiekt, 'pH', 1);
hammerstein.testNonlinearModel(U, a_h, a_pH, b_h, b_pH, obiekt, 'pH', 1);

%% Testy Wiener
U = [repelem((rand(1, obiekt.kk/400) * 30 - 15), 400);
    zeros(1, obiekt.kk);
    repelem((rand(1, obiekt.kk/600) * 30 - 15), 600)];

wiener.testLinearModel(U, a_h, a_Wa4, a_Wb4, b_h, b_Wa4, b_Wb4, obiekt, 'h', 1);
wiener.testNonlinearModel(U, a_h, a_Wa4, a_Wb4, b_h, b_Wa4, b_Wb4, obiekt, 'h', 1);

wiener.testLinearModel(U, a_h, a_Wa4, a_Wb4, b_h, b_Wa4, b_Wb4, obiekt, 'pH', 1);
wiener.testNonlinearModel(U, a_h, a_Wa4, a_Wb4, b_h, b_Wa4, b_Wb4, obiekt, 'pH', 1);

%% DMC(N, Nu, D, D_disturbance, lambda)
S = cell(1, obiekt.dynamic_horizont);

for i = 1:obiekt.dynamic_horizont
    S{i} = [s_h.Q1(i), s_h.Q3(i)
            s_pH.Q1(i), s_pH.Q3(i)];
end

dmc = DMC(100, 10, obiekt.dynamic_horizont, 15);
% Pierwsza iteracja jest wspólna dla wszystkich algorytmów DMC
dmc.dynamic_matrix(S);
dmc.past_matrix(S);

%% Trajektoria zadana
kk = obiekt.kk;

h_levels = zeros(1, ceil(kk/400));
h_levels(1) = 0;
for i = 2:length(h_levels)
    delta = rand()*20 - 10;  % zmiana maks ±5l
    h_levels(i) = max(min(h_levels(i-1) + delta, 10), -10);  % ograniczenie zakresu
end

% --- pH: zmiany schodkowe co 300 próbek ---
pH_levels = zeros(1, ceil(kk/600));
pH_levels(1) = 0;
for i = 2:length(pH_levels)
    delta = rand()*8 - 4;  % zmiana maks ±1
    pH_levels(i) = max(min(pH_levels(i-1) + delta, 4), -4);
end

y_zad(1, :) = repelem(h_levels, 400);
y_zad(2, :) = repelem(pH_levels, 600);

% Ucięcie dokładnie do kk próbek
y_zad(1, :) = y_zad(1, 1:kk);
y_zad(2, :) = y_zad(2, 1:kk);
    
%% Symulacja
[y.analiticH, u.analiticH, E.analiticH, E.u_analiticH, E.y_analiticH] = dmc.dmc_analiticHammerstein(y_zad, a_h, b_h, a_pH, b_pH, obiekt); 
[y.analiticW, u.analiticW, E.analiticW, E.u_analiticW, E.y_analiticW] = dmc.dmc_analiticWiener(y_zad, a_h, b_h, a_Wa4, b_Wa4, a_Wb4, b_Wb4, obiekt); 

t = 0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp;

figure;
plot(t, y_zad(1,:), 'b',...
    t, y.analiticH(1,:), 'r',...
    t, y.analiticW(1,:), ' g'); 
legend('y_zad', 'Hammerstein', 'Wiener');
title('h');
grid on;

figure;
plot(t, y_zad(2, :), 'b',...
    t, y.analiticH(2,:), 'r',...
    t, y.analiticW(2,:), 'g');
legend('y_zad', 'Hammerstein', 'Wiener');
title('pH');
grid on;

%% Wielowątkowy DMC Hammerstein
pool = gcp();

f.analiticHammerstein = parfeval(pool, @dmc.dmc_analiticHammerstein, 5, y_zad, a_h, b_h, a_pH, b_pH, obiekt); 
f.numericHammerstein = parfeval(pool, @dmc.dmc_numericHammerstein, 5, y_zad, a_h, b_h, a_pH, b_pH, obiekt);
f.noHammerstein = parfeval(pool, @dmc.dmc_noHammerstein, 5, y_zad, a_h, b_h, a_pH, b_pH, obiekt);
f.fuzzyHammerstein = parfeval(pool, @dmc.dmc_fuzzyHammerstein, 5, y_zad, a_h, b_h, a_pH, b_pH, obiekt);
f.slH = parfeval(pool, @dmc.dmc_slHammerstein, 5, y_zad, a_h, b_h, a_pH, b_pH, obiekt, hammerstein.linear_fis, 'linear'); 
f.slHammerstein = parfeval(pool, @dmc.dmc_slHammerstein, 5, y_zad, a_h, b_h, a_pH, b_pH, obiekt, hammerstein.nonlinear_fis, 'nonlinear');
f.nplH = parfeval(pool, @dmc.dmc_nplHammerstein, 5, y_zad, a_h, b_h, a_pH, b_pH, obiekt, hammerstein.linear_fis, 'linear');
f.nplHammerstein = parfeval(pool, @dmc.dmc_nplHammerstein, 5, y_zad, a_h, b_h, a_pH, b_pH, obiekt, hammerstein.nonlinear_fis, 'nonlinear');

fprintf("\nFunkcje uruchomione asynchronicznie...\n\n");

[y.analiticHammerstein, u.analiticHammerstein, E.analiticHammerstein, E.u_analiticHammerstein, E.y_analiticHammerstein] = fetchOutputs(f.analiticHammerstein);
% dmc.show_result(y.analiticHammerstein, y_zad, u.analiticHammerstein, E.analiticHammerstein, E.u_analiticHammerstein, E.y_analiticHammerstein, obiekt.kk, 'analitic', 'Hammerstein');

[y.numericHammerstein, u.numericHammerstein, E.numericHammerstein, E.u_numericHammerstein, E.y_numericHammerstein] = fetchOutputs(f.numericHammerstein);
% dmc.show_result(y.numericHammerstein, y_zad, u.numericHammerstein, E.numericHammerstein, E.u_numericHammerstein, E.y_numericHammerstein, obiekt.kk, 'numeric', 'Hammerstein');

[y.noHammerstein, u.noHammerstein, E.noHammerstein, E.u_noHammerstein, E.y_noHammerstein] = fetchOutputs(f.noHammerstein);
% dmc.show_result(y.noHammerstein, y_zad, u.noHammerstein, E.noHammerstein, E.u_noHammerstein, E.y_noHammerstein, obiekt.kk, 'no', 'Hammerstein');

[y.fuzzyHammerstein, u.fuzzyHammerstein, E.fuzzyHammerstein, E.u_fuzzyHammerstein, E.y_fuzzyHammerstein] = fetchOutputs(f.fuzzyHammerstein);
% dmc.show_result(y.fuzzyHammerstein, y_zad, u.fuzzyHammerstein, E.fuzzyHammerstein, E.u_fuzzyHammerstein, E.y_fuzzyHammerstein, obiekt.kk, 'fuzzy', 'Hammerstein');

[y.slHammerstein, u.slHammerstein, E.slHammerstein, E.u_slHammerstein, E.y_slHammerstein] = fetchOutputs(f.slHammerstein);
% dmc.show_result(y.slHammerstein, y_zad, u.slHammerstein, E.slHammerstein, E.u_slHammerstein, E.y_slHammerstein, obiekt.kk, 'SL',' Hammerstein');

[y.nplHammerstein, u.nplHammerstein, E.nplHammerstein, E.u_nplHammerstein, E.y_nplHammerstein] = fetchOutputs(f.nplHammerstein);
% dmc.show_result(y.nplHammerstein, y_zad, u.nplHammerstein, E.nplHammerstein, E.u_nplHammerstein, E.y_nplHammerstein, obiekt.kk, 'NPL', 'Hammerstein');

[y.slH, u.slH, E.slH, E.u_slH, E.y_slH] = fetchOutputs(f.slH);
% dmc.show_result(y.slH, y_zad, u.slH, E.slH, E.u_slH, E.y_slH, obiekt.kk, 'SL', 'Hammerstein');

[y.nplH, u.nplH, E.nplH, E.u_nplH, E.y_nplH] = fetchOutputs(f.nplH);
% dmc.show_result(y.nplH, y_zad, u.nplH, E.nplH, E.u_nplH, E.y_nplH, obiekt.kk, 'NPL', 'Hammerstein');

t = 0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp;

figure;
subplot(2,2,1);
plot(t, y_zad(1,:), 'k-');
hold on;
plot(t, y.analiticHammerstein(1,:), 'b-');
plot(t, y.numericHammerstein(1,:), 'r-');
plot(t, y.slHammerstein(1,:), 'g-');
plot(t, y.nplHammerstein(1,:), 'm-');
plot(t, y.noHammerstein(1,:), 'y-');
plot(t, y.fuzzyHammerstein(1,:), 'c-');
legend('y_{zad}', 'y_{analitic}', ...
    'y_{numeric}', 'y_{slHammerstein}', ...
    'y_{nplHammerstein}', 'y_{no}', 'y_{fuzzy}', 'Location', 'best');
% legend('y_{zad}', 'y_{analitic}', ...
%     'y_{numeric}', 'y_{slHammerstein}', ...
%     'y_{nplHammerstein}', 'y_{fuzzy}', 'Location', 'best');
title(sprintf('Sygnał wyjściowy y(k) - Hammerstein (następniki nieliniowe)'));
xlabel('k');
ylabel('y(k)');
grid on;

subplot(2,2,2);
plot(t, y_zad(1,:), 'k-');
hold on;
plot(t, y.analiticHammerstein(1,:), 'b-');
plot(t, y.numericHammerstein(1,:), 'r-');
plot(t, y.slH(1,:), 'g-');
plot(t, y.nplH(1,:), 'm-');
plot(t, y.noHammerstein(1,:), 'y-');
plot(t, y.fuzzyHammerstein(1,:), 'c-');
legend('y_{zad}', 'y_{analitic}', ...
    'y_{numeric}', 'y_{slHammerstein}', ...
    'y_{nplHammerstein}', 'y_{no}', 'y_{fuzzy}', 'Location', 'best');
% legend('y_{zad}', 'y_{analitic}', ...
%     'y_{numeric}', 'y_{slHammerstein}', ...
%     'y_{nplHammerstein}', 'y_{fuzzy}', 'Location', 'best');
title(sprintf('Sygnał wyjściowy y(k) - Hammerstein (następniki liniowe)'));
xlabel('k');
ylabel('y(k)');
grid on;

subplot(2,2,3);
plot(t, y_zad(2,:), 'k-');
hold on;
plot(t, y.analiticHammerstein(2,:), 'b-');
plot(t, y.numericHammerstein(2,:), 'r-');
plot(t, y.slHammerstein(2,:), 'g-');
plot(t, y.nplHammerstein(2,:), 'm-');
plot(t, y.noHammerstein(2,:), 'y-');
plot(t, y.fuzzyHammerstein(2,:), 'c-');
legend('y_{zad}', 'y_{analitic}', ...
    'y_{numeric}', 'y_{slHammerstein}', ...
    'y_{nplHammerstein}', 'y_{no}', 'y_{fuzzy}', 'Location', 'best');
% legend('y_{zad}', 'y_{analitic}', ...
%     'y_{numeric}', 'y_{slHammerstein}', ...
%     'y_{nplHammerstein}', 'y_{fuzzy}', 'Location', 'best');
title(sprintf('Sygnał wyjściowy y(k) - Hammerstein (następniki nieliniowe)'));
xlabel('k');
ylabel('y(k)');
grid on;

subplot(2,2,4);
plot(t, y_zad(2,:), 'k-');
hold on;
plot(t, y.analiticHammerstein(2,:), 'b-');
plot(t, y.numericHammerstein(2,:), 'r-');
plot(t, y.slH(2,:), 'g-');
plot(t, y.nplH(2,:), 'm-');
plot(t, y.noHammerstein(2,:), 'y-');
plot(t, y.fuzzyHammerstein(2,:), 'c-');
legend('y_{zad}', 'y_{analitic}', ...
    'y_{numeric}', 'y_{slHammerstein}', ...
    'y_{nplHammerstein}', 'y_{no}', 'y_{fuzzy}', 'Location', 'best');
% legend('y_{zad}', 'y_{analitic}', ...
%     'y_{numeric}', 'y_{slHammerstein}', ...
%     'y_{nplHammerstein}', 'y_{fuzzy}', 'Location', 'best');
title(sprintf('Sygnał wyjściowy y(k) - Hammerstein (następniki liniowe)'));
xlabel('k');
ylabel('y(k)');
grid on;

%% Wielowątkowy DMC Wiener
pool = gcp();

f.analiticWiener = parfeval(pool, @dmc.dmc_analiticWiener, 5, y_zad, a_h, b_h, a_Wa4, b_Wa4, a_Wb4, b_Wb4, obiekt); 
f.numericWiener = parfeval(pool, @dmc.dmc_numericWiener, 5, y_zad, a_h, b_h, a_Wa4, b_Wa4, a_Wb4, b_Wb4, obiekt);
f.noWiener = parfeval(pool, @dmc.dmc_noWiener, 5, y_zad, a_h, b_h, a_Wa4, b_Wa4, a_Wb4, b_Wb4, obiekt);
f.fuzzyWiener = parfeval(pool, @dmc.dmc_fuzzyWiener, 5, y_zad, a_h, b_h, a_Wa4, b_Wa4, a_Wb4, b_Wb4, obiekt);
f.slW = parfeval(pool, @dmc.dmc_slWiener, 5, y_zad, a_h, b_h, a_Wa4, b_Wa4, a_Wb4, b_Wb4, obiekt, wiener.linear_fis, 'linear');
f.slWiener = parfeval(pool, @dmc.dmc_slWiener, 5, y_zad, a_h, b_h, a_Wa4, b_Wa4, a_Wb4, b_Wb4, obiekt, wiener.nonlinear_fis, 'nonlinear');
f.nplW = parfeval(pool, @dmc.dmc_nplWiener, 5, y_zad, a_h, b_h, a_Wa4, b_Wa4, a_Wb4, b_Wb4, obiekt, wiener.linear_fis, 'linear');
f.nplWiener = parfeval(pool, @dmc.dmc_nplWiener, 5, y_zad, a_h, b_h, a_Wa4, b_Wa4, a_Wb4, b_Wb4, obiekt, wiener.nonlinear_fis, 'nonlinear');

fprintf("\nFunkcje uruchomione asynchronicznie...\n\n");

[y.analiticWiener, u.analiticWiener, E.analiticWiener, E.u_analiticWiener, E.y_analiticWiener] = fetchOutputs(f.analiticWiener);
% dmc.show_result(y.analiticWiener, y_zad, u.analiticWiener, E.analiticWiener, E.u_analiticWiener, E.y_analiticWiener, obiekt.kk, 'analitic', 'Wiener');

[y.numericWiener, u.numericWiener, E.numericWiener, E.u_numericWiener, E.y_numericWiener] = fetchOutputs(f.numericWiener);
% dmc.show_result(y.numericWiener, y_zad, u.numericWiener, E.numericWiener, E.u_numericWiener, E.y_numericWiener, obiekt.kk, 'numeric', 'Wiener');

[y.noWiener, u.noWiener, E.noWiener, E.u_noWiener, E.y_noWiener] = fetchOutputs(f.noWiener);
% dmc.show_result(y.noWiener, y_zad, u.noWiener, E.noWiener, E.u_noWiener, E.y_noWiener, obiekt.kk, 'no', 'Wiener');

[y.fuzzyWiener, u.fuzzyWiener, E.fuzzyWiener, E.u_fuzzyWiener, E.y_fuzzyWiener] = fetchOutputs(f.fuzzyWiener);
% dmc.show_result(y.fuzzyWiener, y_zad, u.fuzzyWiener, E.fuzzyWiener, E.u_fuzzyWiener, E.y_fuzzyWiener, obiekt.kk, 'fuzzy', 'Wiener');

[y.slWiener, u.slWiener, E.slWiener, E.u_slWiener, E.y_slWiener] = fetchOutputs(f.slWiener);
% dmc.show_result(y.slWiener, y_zad, u.slWiener, E.slWiener, E.u_slWiener, E.y_slWiener, obiekt.kk, 'SL', 'Wiener');

[y.nplWiener, u.nplWiener, E.nplWiener, E.u_nplWiener, E.y_nplWiener] = fetchOutputs(f.nplWiener);
% dmc.show_result(y.nplWiener, y_zad, u.nplWiener, E.nplWiener, E.u_nplWiener, E.y_nplWiener, obiekt.kk, 'NPL', 'Wiener');

[y.slW, u.slW, E.slW, E.u_slW, E.y_slW] = fetchOutputs(f.slW);
% dmc.show_result(y.slW, y_zad, u.slW, E.slW, E.u_slW, E.y_slW, obiekt.kk, 'SL', 'Wiener');

[y.nplW, u.nplW, E.nplW, E.u_nplW, E.y_nplW] = fetchOutputs(f.nplW);
% dmc.show_result(y.nplW, y_zad, u.nplW, E.nplW, E.u_nplW, E.y_nplW, obiekt.kk, 'NPL', 'Wiener');

figure;
subplot(2,2,1);
plot(t, y_zad(1,:), 'k-');
hold on;
plot(t, y.analiticWiener(1,:), 'b-');
plot(t, y.numericWiener(1,:), 'r-');
plot(t, y.slWiener(1,:), 'g-');
plot(t, y.nplWiener(1,:), 'm-');
plot(t, y.noWiener(1,:), 'y-');
plot(t, y.fuzzyWiener(1,:), 'c-');
legend('y_{zad}', 'y_{analitic}', ...
    'y_{numeric}', 'y_{slWiener}', ...
    'y_{nplWiener}', 'y_{no}', 'y_{fuzzy}', 'Location', 'best');
% legend('y_{zad}', 'y_{analitic}', ...
%     'y_{numeric}', 'y_{slWiener}', ...
%     'y_{nplWiener}', 'y_{fuzzy}', 'Location', 'best');
title(sprintf('Sygnał wyjściowy y(k) - Wiener (nastepniki nieliniowe)'));
xlabel('k');
ylabel('y(k)');
grid on;

subplot(2,2,2);
plot(t, y_zad(1,:), 'k-');
hold on;
plot(t, y.analiticWiener(1,:), 'b-');
plot(t, y.numericWiener(1,:), 'r-');
plot(t, y.slW(1,:), 'g-');
plot(t, y.nplW(1,:), 'm-');
plot(t, y.noWiener(1,:), 'y-');
plot(t, y.fuzzyWiener(1,:), 'c-');
legend('y_{zad}', 'y_{analitic}', ...
    'y_{numeric}', 'y_{slWiener}', ...
    'y_{nplWiener}', 'y_{no}', 'y_{fuzzy}', 'Location', 'best');
% legend('y_{zad}', 'y_{analitic}', ...
%     'y_{numeric}', 'y_{slWiener}', ...
%     'y_{nplWiener}', 'y_{fuzzy}', 'Location', 'best');
title(sprintf('Sygnał wyjściowy y(k) - Wiener (nastepniki liniowe)'));
xlabel('k');
ylabel('y(k)');
grid on;

subplot(2,2,3);
plot(t, y_zad(2,:), 'k-');
hold on;
plot(t, y.analiticWiener(2,:), 'b-');
plot(t, y.numericWiener(2,:), 'r-');
plot(t, y.slWiener(2,:), 'g-');
plot(t, y.nplWiener(2,:), 'm-');
plot(t, y.noWiener(2,:), 'y-');
plot(t, y.fuzzyWiener(2,:), 'c-');
legend('y_{zad}', 'y_{analitic}', ...
    'y_{numeric}', 'y_{slWiener}', ...
    'y_{nplWiener}', 'y_{no}', 'y_{fuzzy}', 'Location', 'best');
% legend('y_{zad}', 'y_{analitic}', ...
%     'y_{numeric}', 'y_{slWiener}', ...
%     'y_{nplWiener}', 'y_{fuzzy}', 'Location', 'best');
title(sprintf('Sygnał wyjściowy y(k) - Wiener (nastepniki nieliniowe)'));
xlabel('k');
ylabel('y(k)');
grid on;

subplot(2,2,4);
plot(t, y_zad(2,:), 'k-');
hold on;
plot(t, y.analiticWiener(2,:), 'b-');
plot(t, y.numericWiener(2,:), 'r-');
plot(t, y.slW(2,:), 'g-');
plot(t, y.nplW(2,:), 'm-');
plot(t, y.noWiener(2,:), 'y-');
plot(t, y.fuzzyWiener(2,:), 'c-');
legend('y_{zad}', 'y_{analitic}', ...
    'y_{numeric}', 'y_{slWiener}', ...
    'y_{nplWiener}', 'y_{no}', 'y_{fuzzy}', 'Location', 'best');
% legend('y_{zad}', 'y_{analitic}', ...
%     'y_{numeric}', 'y_{slWiener}', ...
%     'y_{nplWiener}', 'y_{fuzzy}', 'Location', 'best');
title(sprintf('Sygnał wyjściowy y(k) - Wiener (nastepniki liniowe)'));
xlabel('k');
ylabel('y(k)');
grid on;

%% Prezentacja wyników - wysokość h
t = 0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp;

figure;
subplot(2,2,1);
plot(t, y_zad(1,:), 'k-');
hold on;
plot(t, y.analiticHammerstein(1,:), 'b-');
plot(t, y.numericHammerstein(1,:), 'r-');
plot(t, y.slHammerstein(1,:), 'g-');
plot(t, y.nplHammerstein(1,:), 'm-');
plot(t, y.noHammerstein(1,:), 'y-');
plot(t, y.fuzzyHammerstein(1,:), 'c-');
legend('y_{zad}', 'y_{analitic}', ...
    'y_{numeric}', 'y_{slHammerstein}', ...
    'y_{nplHammerstein}', 'y_{no}', 'y_{fuzzy}', 'Location', 'best');
title(sprintf('Sygnał wyjściowy y(k) - Hammerstein'));
xlabel('k');
ylabel('y(k)');
grid on;

subplot(2,2,2);
plot(t, y_zad(1,:), 'k-');
hold on;
plot(t, y.analiticWiener(1,:), 'b-');
plot(t, y.numericWiener(1,:), 'r-');
plot(t, y.slWiener(1,:), 'g-');
plot(t, y.nplWiener(1,:), 'm-');
plot(t, y.noWiener(1,:), 'y-');
plot(t, y.fuzzyWiener(1,:), 'c-');
legend('y_{zad}', 'y_{analitic}', ...
    'y_{numeric}', 'y_{slWiener}', ...
    'y_{nplWiener}', 'y_{no}', 'y_{fuzzy}', 'Location', 'best');
title(sprintf('Sygnał wyjściowy y(k) - Wiener'));
xlabel('k');
ylabel('y(k)');
grid on;

subplot(2,2,3);
plot(t, y_zad(1,:), 'k-');
hold on;
plot(t, y.analiticHammerstein(1,:), 'b-');
plot(t, y.numericHammerstein(1,:), 'r-');
plot(t, y.slH(1,:), 'g-');
plot(t, y.nplH(1,:), 'm-');
plot(t, y.noHammerstein(1,:), 'y-');
plot(t, y.fuzzyHammerstein(1,:), 'c-');
legend('y_{zad}', 'y_{analitic}', ...
    'y_{numeric}', 'y_{slHammerstein}', ...
    'y_{nplHammerstein}', 'y_{no}', 'y_{fuzzy}', 'Location', 'best');
title(sprintf('Sygnał wyjściowy y(k) - Hammerstein'));
xlabel('k');
ylabel('y(k)');
grid on;

subplot(2,2,4);
plot(t, y_zad(1,:), 'k-');
hold on;
plot(t, y.analiticWiener(1,:), 'b-');
plot(t, y.numericWiener(1,:), 'r-');
plot(t, y.slW(1,:), 'g-');
plot(t, y.nplW(1,:), 'm-');
plot(t, y.noWiener(1,:), 'y-');
plot(t, y.fuzzyWiener(1,:), 'c-');
legend('y_{zad}', 'y_{analitic}', ...
    'y_{numeric}', 'y_{slWiener}', ...
    'y_{nplWiener}', 'y_{no}', 'y_{fuzzy}', 'Location', 'best');
title(sprintf('Sygnał wyjściowy y(k) - Wiener'));
xlabel('k');
ylabel('y(k)');
grid on;

%% Prezentacja wyników - stężenie pH
figure;
subplot(2,2,1);
plot(t, y_zad(2,:), 'k-');
hold on;
plot(t, y.analiticHammerstein(2,:), 'b-');
plot(t, y.numericHammerstein(2,:), 'r-');
plot(t, y.slHammerstein(2,:), 'g-');
plot(t, y.nplHammerstein(2,:), 'm-');
plot(t, y.noHammerstein(2,:), 'y-');
plot(t, y.fuzzyHammerstein(2,:), 'c-');
legend('y_{zad}', 'y_{analitic}', ...
    'y_{numeric}', 'y_{slHammerstein}', ...
    'y_{nplHammerstein}', 'y_{no}', 'y_{fuzzy}', 'Location', 'best');
title(sprintf('Sygnał wyjściowy y(k) - Hammerstein'));
xlabel('k');
ylabel('y(k)');
grid on;

subplot(2,2,2);
plot(t, y_zad(2,:), 'k-');
hold on;
plot(t, y.analiticWiener(2,:), 'b-');
plot(t, y.numericWiener(2,:), 'r-');
plot(t, y.slWiener(2,:), 'g-');
plot(t, y.nplWiener(2,:), 'm-');
plot(t, y.noWiener(2,:), 'y-');
plot(t, y.fuzzyWiener(2,:), 'c-');
legend('y_{zad}', 'y_{analitic}', ...
    'y_{numeric}', 'y_{slWiener}', ...
    'y_{nplWiener}', 'y_{no}', 'y_{fuzzy}', 'Location', 'best');
title(sprintf('Sygnał wyjściowy y(k) - Wiener'));
xlabel('k');
ylabel('y(k)');
grid on;

subplot(2,2,3);
plot(t, y_zad(2,:), 'k-');
hold on;
plot(t, y.analiticHammerstein(2,:), 'b-');
plot(t, y.numericHammerstein(2,:), 'r-');
plot(t, y.slH(2,:), 'g-');
plot(t, y.nplH(2,:), 'm-');
plot(t, y.noHammerstein(2,:), 'y-');
plot(t, y.fuzzyHammerstein(2,:), 'c-');
legend('y_{zad}', 'y_{analitic}', ...
    'y_{numeric}', 'y_{slHammerstein}', ...
    'y_{nplHammerstein}', 'y_{no}', 'y_{fuzzy}', 'Location', 'best');
title(sprintf('Sygnał wyjściowy y(k) - Hammerstein'));
xlabel('k');
ylabel('y(k)');
grid on;

subplot(2,2,4);
plot(t, y_zad(2,:), 'k-');
hold on;
plot(t, y.analiticWiener(2,:), 'b-');
plot(t, y.numericWiener(2,:), 'r-');
plot(t, y.slW(2,:), 'g-');
plot(t, y.nplW(2,:), 'm-');
plot(t, y.noWiener(2,:), 'y-');
plot(t, y.fuzzyWiener(2,:), 'c-');
legend('y_{zad}', 'y_{analitic}', ...
    'y_{numeric}', 'y_{slWiener}', ...
    'y_{nplWiener}', 'y_{no}', 'y_{fuzzy}', 'Location', 'best');
title(sprintf('Sygnał wyjściowy y(k) - Wiener'));
xlabel('k');
ylabel('y(k)');
grid on;

% delete(pool);