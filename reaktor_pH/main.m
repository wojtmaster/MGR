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
[a_h, b_h, s_h] = obiekt.mse('h');
[a_pH, b_pH, s_pH] = obiekt.mse('pH');

S = cell(1, obiekt.dynamic_horizont);
S_disturbance = cell(1, obiekt.dynamic_horizont);

for i = 1:obiekt.dynamic_horizont
    S{i} = [s_h.Q1(1,i), s_h.Q3(1,i)
            s_pH.Q1(2,i), s_pH.Q3(2,i)];
    S_disturbance{i} = [s_h.Q2(1,i); s_pH.Q2(2,i)];
end

%% Hammerstein 
hammerstein = Hammerstein(@(u, kk) obiekt.modifiedEuler(u, kk));
hammerstein.linearFuzzy();
hammerstein.nonlinearFuzzy();

%% Wiener
wiener = Wiener(a_h, a_pH, b_h, b_pH, obiekt);
wiener.linearFuzzy();
wiener.nonlinearFuzzy();

%% Check fuzzy static
U_min = -15;
U_max = 15;
U = [linspace(U_min, U_max, 100);
    linspace(U_min, U_max, 100)];
[Q1_grid, Q3_grid] = meshgrid(U(1,:), U(2,:));

% % Dane wejściowe
% X = [Q1_grid(:), Q3_grid(:)]';
% Y = wiener.pH(:)';
% 
% % Sieć neuronowa (feedforward)
% net = fitnet([10 10], 'trainlm');  % 2 warstwy ukryte po 10 neuronów
% net.trainParam.showWindow = true;
% net = train(net, X, Y);
% 
% % Predykcja
% Y_pred = net(X)';
% pH_pred = reshape(Y_pred, size(Q1_grid));
% 
% % Rysowanie
% figure;
% surf(Q1_grid, Q3_grid, pH_pred);
% xlabel('Q_1'); ylabel('Q_3'); zlabel('pH');
% title('Charakterystyka statyczna - sieć neuronowa');
% shading interp; colormap parula;

Y_out = zeros(100,100);
for i = 1:100
    disp(i);
    for j = 1:100
        Y_out(i,j) = evalfis(hammerstein.linear_fis.pH, [U(1,j) U(2,i)]);
        % Y_out(i,j) = evalfis(wiener.linear_fis.pH, wiener.Y_pH(i,j));
        % Y_out(i,j) = evalfis(hammerstein.nonlinear_fis.pH, tanh((U(1,j)-U(2,i))/15));
        % Y_out(i,j) = evalfis(wiener.nonlinear_fis.pH, tanh(wiener.Y_pH(i,j)/30));
    end
end

figure;
surf(Q1_grid, Q3_grid, Y_out);
xlabel('Q_1 [ml/s]');
ylabel('Q_3 [ml/s]');
zlabel('pH');
title(sprintf('Wpływ dopływów Q_1 oraz Q_3 na stężenie substancji pH\nHammerstein- następniki liniowe'));
shading interp;
colorbar;
view(-45, 30);

%% Testy Hammerstein
U = [repelem((rand(1, obiekt.kk/300) * 30 - 15), 300);
    zeros(1, obiekt.kk);
    repelem((rand(1, obiekt.kk/700) * 30 - 15), 700)];

% hammerstein.testLinearModel(U, a_h, a_pH, b_h, b_pH, obiekt, 'h', 1);
% hammerstein.testNonlinearModel(U, a_h, a_pH, b_h, b_pH, obiekt, 'h', 1);

% hammerstein.testLinearModel(U, a_h, a_pH, b_h, b_pH, obiekt, 'pH', 1);
hammerstein.testNonlinearModel(U, a_h, a_pH, b_h, b_pH, obiekt, 'pH', 1);

%% Testy Wiener
% U = [repelem((rand(1, obiekt.kk/300) * 30 - 15), 300);
%     zeros(1, obiekt.kk);
%     repelem((rand(1, obiekt.kk/700) * 30 - 15), 700)];

% wiener.testLinearModel(U, a_h, a_pH, b_h, b_pH, obiekt, 'h', 1);
% wiener.testNonlinearModel(U, a_h, a_pH, b_h, b_pH, obiekt, 'h', 1);

wiener.testLinearModel(U, a_h, a_pH, b_h, b_pH, obiekt, 'pH', 1);
wiener.testNonlinearModel(U, a_h, a_pH, b_h, b_pH, obiekt, 'pH', 1);

%% DMC(N, Nu, D, D_disturbance, lambda)
dmc = DMC(100, 4, 120, 120, 1);
% Pierwsza iteracja jest wspólna dla wszystkich algorytmów DMC
dmc.dynamic_matrix(S);
dmc.past_matrix(S);
dmc.matrix_disturbance(S_disturbance);

%% Test
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