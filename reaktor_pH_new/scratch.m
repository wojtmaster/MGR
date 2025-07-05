%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obiekt = Obiekt();
obiekt.linearization(16.6, 0.55, 15.6);
U = [linspace(1.6, 31.6, 100)
    linspace(0.6, 30.6, 100)];
[Q1_grid, Q3_grid] = meshgrid(U(1,:), U(2,:));
Q2 = 0.55;
h = ((Q1_grid + Q2 + Q3_grid) / obiekt.C_V).^2 - obiekt.h_0;
Wa4 = (obiekt.W_a1*Q1_grid + obiekt.W_a2*Q2 + obiekt.W_a3*Q3_grid)./(Q1_grid+Q2+Q3_grid);
Wb4 = (obiekt.W_b1*Q1_grid + obiekt.W_b2*Q2 + obiekt.W_b3*Q3_grid)./(Q1_grid+Q2+Q3_grid);

for i = 1:100
    for j = 1:100
        pH(i,j) = obiekt.pH_calc(Wa4(i,j), Wb4(i,j)) - obiekt.pH_0;
    end
end

figure;
surf(Q1_grid, Q3_grid, h);

figure;
surf(Q1_grid, Q3_grid, Wa4);

figure;
surf(Q1_grid, Q3_grid, Wb4);

figure;
surf(Q1_grid, Q3_grid, pH);

%% NONLINEAR HAMMERSTEIN MODEL
clear all;
obiekt = Obiekt();
[a_1, a_2, a_3, b_1, b_2, b_3, G_z] = obiekt.linearization(16.6, 0.55, 15.6);

% Zakres sterowania
U_min = -15;
U_max = 15;
U_center = linspace(U_min, U_max, 3); % Środek zbiorów

for i = 1:50
    U_tmp = [repelem((rand(1, obiekt.kk/400) * 30 - 15), 400);
             zeros(1, obiekt.kk);
             repelem((rand(1, obiekt.kk/400) * 30 - 15), 400)];
    [Y_tmp, ~] = obiekt.modifiedEuler(U_tmp, obiekt.kk);

    U_Q1(i, :) = U_tmp(1, :);
    U_Q3(i, :) = U_tmp(3, :);
    Y_train(i, :) = Y_tmp(1, :);
end
t = (0:length(U_Q1(1,:))-1) * obiekt.Tp;

%% FMINSEARCH
close all;
fis = sugfis('Name', 'Hammerstein', 'Type', 'sugeno');
rules_number = 6;

fis = addInput(fis, [U_min U_max], 'Name', 'u1');
fis = addInput(fis, [U_min U_max], 'Name', 'u2');

% Definiowanie funkcji przynależności (gaussmf)
for i = 1:length(U_center) 
    fis = addMF(fis, 'u1', 'gaussmf', [6, U_center(i)]);
    fis = addMF(fis, 'u2', 'gaussmf', [6, U_center(i)]);
end

% Definiowanie wyjścia i początkowych następników (a_i * u + b_i)
fis = addOutput(fis, [U_min U_max], 'Name', 'u_fuzzy');

% Początkowe współczynniki (a_i, b_i, c_i)
a_param = ones(1,rules_number)*0.01;
b_param = ones(1,rules_number)*15;
c_param = ones(1,rules_number)*15;
% a_param = [0.2773    0.9531  -0.0001    0.9504    0.3924    0.0766];
% b_param = [3.8845    4.7977    1.4201    3.7108    3.5950    0.0135];
% c_param = [3.5019    5.0929    3.1197    4.6958    3.5080  -20.5576];

% Dodanie reguł TS w postaci liniowej
for i = 1:rules_number
    fis = addMF(fis, 'u_fuzzy', 'linear', [a_param(i), b_param(i), c_param(i)]);
end

% Reguły Takagi-Sugeno: [inputMF, outputMF, weight]
ruleList = [1 1 1 1 1;  
            1 2 2 1 1;
            2 1 2 1 1;
            1 3 3 1 1;
            3 1 3 1 1;
            2 2 4 1 1;
            2 3 5 1 1;
            3 2 5 1 1;
            3 3 6 1 1];

% Dodanie reguł do systemu
fis = addRule(fis, ruleList);

% Optymalizacja przy pomocy fminsearch
initial_params = [a_param, b_param c_param];  % Początkowe wartości a_param i b_param
options = optimset('Display', 'iter', 'MaxFunEvals', 2500, 'MaxIter', 2500); % Opcje optymalizacji
optimal_params = fminsearch(@(params) nonlinearCoeff(params, fis, U_Q1, U_Q3, Y_train, U_center, rules_number, a_1, b_1.Q_1, b_1.Q_3), initial_params, options);

% Po optymalizacji
a_optimal = optimal_params(1:rules_number);
b_optimal = optimal_params(rules_number+1:2*rules_number);
c_optimal = optimal_params(2*rules_number+1:end);

% Wyświetlanie wyników optymalizacji
fprintf('Optymalne parametry a: \n');
disp(a_optimal);
fprintf('Optymalne parametry b: \n');
disp(b_optimal);
fprintf('Optymalne parametry c: \n');
disp(c_optimal);

%% 2. Tworzenie początkowego systemu rozmytego TS
% a_optimal = [0.0035    0.0013    0.0211    0.0469   -0.0046    3.1275];
% b_optimal = [134.4978  473.1631  178.4632   94.5173  10.0297   25.2278];
% c_optimal = [130.9950  165.6131  327.5411   89.8189   48.8507 -15.9479];
for i = 1:rules_number
    fis.Outputs.MembershipFunctions(i).Parameters(1) = a_optimal(i);
    fis.Outputs.MembershipFunctions(i).Parameters(2) = b_optimal(i);
    fis.Outputs.MembershipFunctions(i).Parameters(3) = c_optimal(i);
end

% Symulacja modelu Hammersteina
U_fuzzy = zeros(1, obiekt.kk); % Przepuszczenie przez model TS
Y_out = zeros(1, obiekt.kk);
rule_to_param = [1 2 2 3 3 1 5 5 4];

U = [repelem((rand(1, obiekt.kk/400) * 30 - 15), 400);
    zeros(1, obiekt.kk);
    repelem((rand(1, obiekt.kk/400) * 30 - 15), 400)];
[Y_real, Y_lin] = obiekt.modifiedEuler(U, obiekt.kk); % Symulacja rzeczywistego układu

for k = 2:obiekt.kk
    [~, degrees] = evalfis(fis, [U(1,k), U(3,k)]);
    output = 0;
    w = 0;
    for i = 1:length(fis.Rules)
        index = rule_to_param(i);
        output = output + (degrees(i,1)*degrees(i,2))*fis.Outputs.MembershipFunctions(index).Parameters(1)*...
            (sinh(U(1,k)/fis.Outputs.MembershipFunctions(index).Parameters(2)) + sinh(U(3,k)/fis.Outputs.MembershipFunctions(index).Parameters(3)));
        w = w + (degrees(i,1)*degrees(i,2));
    end
    U_fuzzy(k) = output / w;
    Y_out(k) = - a_1*Y_out(k-1) + b_1.Q_1 * U_fuzzy(k-1) + b_1.Q_2 * U(2, k-1) + b_1.Q_3 * U_fuzzy(k-1);
end

% Wizualizacja wyników
figure;
plot(t, Y_real(1,:), 'b', t, Y_lin(1,:), 'g', t, Y_out, 'r');
legend('Euler', 'Euler liniowy', 'Hammerstein (optymalny TS)', 'Location', 'southwest');
title('Porównanie wyjścia układu rzeczywistego i modelu');
grid on;

E_lin = sum((Y_real(1,:) - Y_lin(1,:)).^2);
E_out = sum((Y_real(1,:) - Y_out).^2);
fprintf("\nE_lin = %.3f\n", E_lin);
fprintf("E_out = %.3f\n", E_out);

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
Y_fuzzy_wa4 = evalfis(hammerstein.linear_fis.Wa4, [U(1,:)' U(3,:)']);
Y_fuzzy_wb4 = evalfis(hammerstein.linear_fis.Wb4, [U(1,:)' U(3,:)']);

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

%% Do obliczania macierzy K w DMC 
% [nRows, nCols] = size(obj.K);
% 
% % Sprawdzenie, czy da się podzielić na bloki 2x2
% if mod(nRows,2) ~= 0 || mod(nCols,2) ~= 0
%     error('Wymiary macierzy muszą być podzielne przez 2.');
% end
% 
% % Wektory rozmiarów bloków
% rowBlocks = repmat(2, 1, nRows/2);  % np. [2 2 2 2]
% colBlocks = repmat(2, 1, nCols/2);  % np. [2 2 ... x100]
% 
% % Podział na komórki 2x2
% K_cell = mat2cell(obj.K, rowBlocks, colBlocks);