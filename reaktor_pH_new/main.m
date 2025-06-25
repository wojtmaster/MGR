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

[a_h, b_h, s_h] = obiekt.fopdtModel('h', 1);

[a_pH, b_pH, s_pH] = obiekt.tfestModel();

[a_Wa4, b_Wa4, s_Wa4] = obiekt.fopdtModel('W_{a4}', 3);
[a_Wb4, b_Wb4, s_Wb4] = obiekt.fopdtModel('W_{b4}', 4);

%% Charakterystyki statyczne Wa4 i Wb4
U = [linspace(1.6, 31.6, 100)
    linspace(0.6, 30.6, 100)];
[Q1_grid, Q3_grid] = meshgrid(U(1,:), U(2,:));
Wa4 = (obiekt.W_a1*Q1_grid + obiekt.W_a2*obiekt.Q_20 + obiekt.W_a3*Q3_grid)./(Q1_grid+obiekt.Q_20+Q3_grid);
Wb4 = (obiekt.W_b1*Q1_grid + obiekt.W_b2*obiekt.Q_20 + obiekt.W_b3*Q3_grid)./(Q1_grid+obiekt.Q_20+Q3_grid);

for i = 1:100
    for j = 1:100
        pH(i,j) = obiekt.pH_calc(Wa4(i,j), Wb4(i,j)) - obiekt.pH_0;
    end
end

figure;
surf(Q1_grid, Q3_grid, Wa4);
xlabel('Q_1 [ml/s]');
ylabel('Q_3 [ml/s]');
zlabel('Wa4');
title('Wpływ dopływów Q_1 oraz Q_3 na wysokość słupa cieczy Wa4');
shading interp;
colorbar;

figure;
surf(Q1_grid, Q3_grid, Wb4);
xlabel('Q_1 [ml/s]');
ylabel('Q_3 [ml/s]');
zlabel('Wb4');
title('Wpływ dopływów Q_1 oraz Q_3 na wysokość słupa cieczy Wb4');
shading interp;
colorbar;

Wa4 = Wa4 - obiekt.W_a40;
Wb4 = Wb4 - obiekt.W_b40;

figure;
surf(Wb4, Wa4, pH);
xlabel('Q_1 [ml/s]');
ylabel('Q_3 [ml/s]');
zlabel('Wb4');
title('Wpływ dopływów Q_1 oraz Q_3 na wysokość słupa cieczy Wb4');
shading interp;
colorbar;

U_min = -15;
U_max = 15;
U = [linspace(U_min, U_max, 100);
    linspace(U_min, U_max, 100)];

[Q1_grid, Q3_grid] = meshgrid(U(1,:), U(2,:));  % kombinacje

% Dane wejściowe
X1 = Q1_grid';
X2 = Q3_grid';
X_Wa4 = [X1(:) X2(:)];
X_Wb4 = [X1(:) X2(:)];

Y_Wa4 = Wa4';
Y_Wa4 = Y_Wa4(:);

Y_Wb4 = Wb4';
Y_Wb4 = Y_Wb4(:);

% Generowanie systemu rozmytego
fis_Wa4 = genfis1([X_Wa4 Y_Wa4], 5, 'gaussmf', 'linear');
fis_Wb4 = genfis1([X_Wb4 Y_Wb4], 5, 'gaussmf', 'linear');

% Trening ANFIS
options_1 = anfisOptions('InitialFIS', fis_Wa4, 'EpochNumber', 100, ...
    'DisplayANFISInformation', 0, 'DisplayErrorValues', 0);
options_2 = anfisOptions('InitialFIS', fis_Wb4, 'EpochNumber', 100, ...
    'DisplayANFISInformation', 0, 'DisplayErrorValues', 0);
fis_trained_Wa4 = anfis([X_Wa4 Y_Wa4], options_1);
fis_trained_Wb4 = anfis([X_Wb4 Y_Wb4], options_2);

% Y_Wa = zeros(100,100);
% Y_Wb = zeros(100,100);
% Y_pH = zeros(100, 100);
% for i = 1:100
%     for j = 1:100
%         Y_Wa(i,j) = evalfis(fis_trained_Wa4, [U(1,j) U(2,i)]);
%         Y_Wb(i,j) = evalfis(fis_trained_Wb4, [U(1,j) U(2,i)]);
%         Y_pH(i,j) = obiekt.pH_calc(Y_Wa(i,j)+obiekt.W_a40, Y_Wb(i,j)+obiekt.W_b40) - obiekt.pH_0;
%     end
% end
% 
% % Rysunek
% figure;
% surf(Q1_grid, Q3_grid, Y_pH);
% xlabel('Q_1'); ylabel('Q_3'); zlabel('pH');
% title('Dokładniejsza charakterystyka statyczna TS (genfis1 + ANFIS)');
% shading interp; colormap turbo;

%% Hammerstein 
hammerstein = Hammerstein(obiekt);
hammerstein.linearFuzzy();
hammerstein.nonlinearFuzzy();

%% Wiener
wiener = Wiener(a_h, a_Wa4, a_Wb4, b_h, b_Wa4, b_Wb4, obiekt);
wiener.linearFuzzy();
wiener.nonlinearFuzzy();

%% Symulacja
U = [repelem((rand(1, obiekt.kk/300) * 30 - 15), 300);
    zeros(1, obiekt.kk);
    repelem((rand(1, obiekt.kk/700) * 30 - 15), 700)];

[Y_real, ~] = obiekt.modifiedEuler(U, obiekt.kk);

y_Wa4.Q1 = zeros(1, obiekt.kk);
y_Wa4.Q3 = zeros(1, obiekt.kk);
y_Wb4.Q1 = zeros(1, obiekt.kk);
y_Wb4.Q3 = zeros(1, obiekt.kk);
Wa4 = zeros(1, obiekt.kk);
Wb4 = zeros(1, obiekt.kk);

for k = 3:obiekt.kk
    % y_Wa4.Q1(k) = - a_Wa4*y_Wa4.Q1(k-1:-1:k-2)' + b_Wa4*U(1, k-1:-1:k-2)';
    % y_Wa4.Q3(k) = - a_Wa4*y_Wa4.Q3(k-1:-1:k-2)' - b_Wa4*U(3, k-1:-1:k-2)';
    % 
    % y_Wb4.Q1(k) = - a_Wb4*y_Wb4.Q1(k-1:-1:k-2)' + b_Wb4*U(1, k-1:-1:k-2)';
    % y_Wb4.Q3(k) = - a_Wb4*y_Wb4.Q3(k-1:-1:k-2)' + b_Wb4*U(3, k-1:-1:k-2)';

    y_Wa4.Q1(k) = - a_Wa4*y_Wa4.Q1(k-1)' + b_Wa4*U(1, k-1)';
    y_Wa4.Q3(k) = - a_Wa4*y_Wa4.Q3(k-1)' - b_Wa4*U(3, k-1)';

    y_Wb4.Q1(k) = - a_Wb4*y_Wb4.Q1(k-1)' + b_Wb4*U(1, k-1)';
    y_Wb4.Q3(k) = - a_Wb4*y_Wb4.Q3(k-1)' + b_Wb4*U(3, k-1)';
    
    Wa4(k) = (obiekt.W_a1*(obiekt.Q_10+y_Wa4.Q1(k)) + obiekt.W_a2*obiekt.Q_20 + obiekt.W_a3*(obiekt.Q_30 - y_Wa4.Q3(k)))./(obiekt.Q_10+obiekt.Q_20+obiekt.Q_30+y_Wa4.Q1(k)-y_Wa4.Q3(k)) - obiekt.W_a40;
    Wb4(k) = (obiekt.W_b1*(obiekt.Q_10-y_Wb4.Q1(k)) + obiekt.W_b2*obiekt.Q_20 + obiekt.W_b3*(obiekt.Q_30- y_Wb4.Q3(k)))./(obiekt.Q_10+obiekt.Q_20+obiekt.Q_30-y_Wb4.Q1(k)- y_Wb4.Q3(k)) - obiekt.W_b40;
    pH(k) = obiekt.pH_calc(Wa4(k) + obiekt.W_a40, Wb4(k) + obiekt.W_b40) - obiekt.pH_0;
end

figure;
plot(Y_real(3,:), 'b');
hold on;
plot(Wa4, 'r');
title('Wa4');

figure;
plot(Y_real(4,:), 'b');
hold on;
plot(Wb4, 'r');
title('Wb4');

figure;
plot(Y_real(2,:), 'b');
hold on;
plot(pH, 'r');
title('pH');

%% Check fuzzy static
U_min = -15;
U_max = 15;
U = [linspace(U_min, U_max, 100);
    linspace(U_min, U_max, 100)];
[Q1_grid, Q3_grid] = meshgrid(U(1,:), U(2,:));

Y_Hout = zeros(100,100);
Y_Aout = zeros(100,100);
Y_Bout = zeros(100,100);
for i = 1:100
    disp(i);
    for j = 1:100
        % Y_out(i,j) = evalfis(hammerstein.linear_fis.h, [U(1,j) U(2,i)]);
        % Y_out(i,j) = evalfis(hammerstein.nonlinear_fis.h, [sinh(U(1,j)/15) sinh(U(2,i)/15)]);
        % Y_out(i,j) = evalfis(wiener.linear_fis.h, wiener.Y_h(i,j));
        % Y_out(i,j) = evalfis(wiener.nonlinear_fis.h, sinh(wiener.Y_h(i,j)/30));

        % Y_Hout(i,j) = evalfis(wiener.linear_fis.h, [wiener.Y_h.Q1(i,j), wiener.Y_h.Q3(i,j)]);
        % Y_Aout(i,j) = evalfis(wiener.linear_fis.Wa4, [wiener.Y_Wa4.Q1(i,j), wiener.Y_Wa4.Q3(i,j)]);
        % Y_Bout(i,j) = evalfis(wiener.linear_fis.Wb4, [wiener.Y_Wb4.Q1(i,j), wiener.Y_Wb4.Q3(i,j)]);

        Y_Hout(i,j) = evalfis(wiener.nonlinear_fis.h, [sinh(wiener.Y_h.Q1(i,j)/15), sinh(wiener.Y_h.Q3(i,j)/15)]);
        Y_Aout(i,j) = evalfis(wiener.nonlinear_fis.Wa4, [tanh(wiener.Y_Wa4.Q1(i,j)/15), tanh(wiener.Y_Wa4.Q3(i,j)/15)]);
        Y_Bout(i,j) = evalfis(wiener.nonlinear_fis.Wb4, [sinh(wiener.Y_Wb4.Q1(i,j)/15), sinh(wiener.Y_Wb4.Q3(i,j)/15)]);
    end
end

figure;
surf(Q1_grid, Q3_grid, Y_Hout);
xlabel('Q_1 [ml/s]');
ylabel('Q_3 [ml/s]');
zlabel('h');
title(sprintf('h'));
shading interp;
colorbar;

figure;
surf(Q1_grid, Q3_grid, Y_Aout);
xlabel('Q_1 [ml/s]');
ylabel('Q_3 [ml/s]');
zlabel('pH');
title(sprintf('Wa4'));
shading interp;
colorbar;

figure;
surf(Q1_grid, Q3_grid, Y_Bout);
xlabel('Q_1 [ml/s]');
ylabel('Q_3 [ml/s]');
zlabel('pH');
title(sprintf('Wb4'));
shading interp;
colorbar;

%% Losowa sekwencja
U = [repelem((rand(1, obiekt.kk/300) * 30 - 15), 300);
    zeros(1, obiekt.kk);
    repelem((rand(1, obiekt.kk/700) * 30 - 15), 700)];

[Y_real, Y_lin] = obiekt.modifiedEuler(U, obiekt.kk);
t = 0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp; 
v_wa4 = zeros(1, obiekt.kk);
v_wa4(1:2) = Y_real(3,1:2);
v_wb4 = zeros(1, obiekt.kk);
v_wb4(1:2) = Y_real(4,1:2);
du = 0.1;
Y_fuzzy_wa4 = evalfis(fis_trained_Wa4, [U(1,:)' U(3,:)']);
Y_fuzzy_wb4 = evalfis(fis_trained_Wb4, [U(1,:)' U(3,:)']);

% Y_fuzzy_wa4_q1 = evalfis(fis_trained_Wa4, [U(1,:)'-du U(3,:)']);
% Y_fuzzy_wb4_q1 = evalfis(fis_trained_Wb4, [U(1,:)'-du U(3,:)']);

% Y_fuzzy_wa4_q3 = evalfis(fis_trained_Wa4, [U(1,:)' U(3,:)'-du]);
% Y_fuzzy_wb4_q3 = evalfis(fis_trained_Wb4, [U(1,:)' U(3,:)'-du]);

Y_out = zeros(1, obiekt.kk);
Y_out(1:2) = Y_real(2, 1:2);

for k = 3:obiekt.kk
    % dydq1_wa4 = (Y_fuzzy_wa4(k-1) - Y_fuzzy_wa4_q1(k-1)) / du;
    % dydq1_wb4 = (Y_fuzzy_wb4(k-1) - Y_fuzzy_wb4_q1(k-1)) / du;
    % dydq3_wa4 = (Y_fuzzy_wa4(k-1) - Y_fuzzy_wa4_q3(k-1)) / du;
    % dydq3_wb4 = (Y_fuzzy_wb4(k-1) - Y_fuzzy_wb4_q3(k-1)) / du;
    
    K_wa4 = Y_fuzzy_wa4(k-1) / (U(1, k-1) - U(3, k-1));
    K_wb4 = Y_fuzzy_wb4(k-1) / (- U(1, k-1) - U(3, k-1));
    
    % disp(dydq1_wa4);
    % disp(dydq1_wb4);
    % disp(dydq3_wa4);
    % disp(dydq3_wb4);

    v_wa4(k) = -a_Wa4 * v_wa4(k-1) + K_wa4*b_Wa4 * U(1, k-1) - K_wa4*b_Wa4 * U(3, k-1);
    v_wb4(k) = -a_Wb4 * v_wb4(k-1) + K_wb4*b_Wb4 * U(1, k-1) + K_wb4*b_Wb4 * U(3, k-1);
    Y_out(k) = obiekt.pH_calc(obiekt.W_a40 + v_wa4(k), obiekt.W_b40 + v_wb4(k)) - obiekt.pH_0;
end

figure;
plot(t, Y_real(2, :), 'b', t, Y_lin(2, :), 'g', t, Y_out, 'r', 'LineWidth', 1.5);
% plot(t, Y_real(3, :), 'b', t, Y_lin(3, :), 'g', 'LineWidth', 1.5);

%% Testy Hammerstein
U = [repelem((rand(1, obiekt.kk/300) * 30 - 15), 300);
    zeros(1, obiekt.kk);
    zeros(1, obiekt.kk)];

hammerstein.testLinearModel(U, a_h, a_pH, b_h, b_pH, obiekt, 'h', 1);
hammerstein.testNonlinearModel(U, a_h, a_pH, b_h, b_pH, obiekt, 'h', 1);

hammerstein.testLinearModel(U, a_h, a_pH, b_h, b_pH, obiekt, 'pH', 1);
hammerstein.testNonlinearModel(U, a_h, a_pH, b_h, b_pH, obiekt, 'pH', 1);

%% Testy Wiener
U = [repelem((rand(1, obiekt.kk/300) * 30 - 15), 300);
    zeros(1, obiekt.kk);
    repelem((rand(1, obiekt.kk/700) * 30 - 15), 700)];

wiener.testLinearModel(U, a_h, a_Wa4, a_Wb4, b_h, b_Wa4, b_Wb4, obiekt, 'h', 1);
wiener.testNonlinearModel(U, a_h, a_Wa4, a_Wb4, b_h, b_Wa4, b_Wb4, obiekt, 'h', 1);

wiener.testLinearModel(U, a_h, a_Wa4, a_Wb4, b_h, b_Wa4, b_Wb4, obiekt, 'pH', 1);
wiener.testNonlinearModel(U, a_h, a_Wa4, a_Wb4, b_h, b_Wa4, b_Wb4, obiekt, 'pH', 1);

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