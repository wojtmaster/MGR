%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




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

%% 
%% 2. Tworzenie początkowego systemu rozmytego TS
close all;
fis_2 = sugfis('Name', 'F1_Hammerstein', 'Type', 'sugeno');
fis_2 = addInput(fis_2, [U_min U_max], 'Name', 'U');

% Definiowanie funkcji przynależności (gaussmf)
fis_2 = addMF(fis_2, 'U', 'gaussmf', [20, U_center(1)]);
fis_2 = addMF(fis_2, 'U', 'gaussmf', [20, U_center(2)]);
fis_2 = addMF(fis_2, 'U', 'gaussmf', [20, U_center(3)]);

% Definiowanie wyjścia i początkowych następników (a_i * u + b_i)
fis_2 = addOutput(fis_2, [U_min U_max], 'Name', 'F1');

% Początkowe współczynniki (a_i, b_i)
a_param = [76.4472   72.8656  100.0000]; % Można dobrać inaczej
b_param = [100.0000   71.1456   88.1261];

% Dodanie reguł TS w postaci liniowej
for i = 1:length(U_center)
    fis_2 = addMF(fis_2, 'F1', 'linear', [a_param(i), b_param(i)]);
end

% Reguły Takagi-Sugeno: [inputMF, outputMF, weight]
ruleList = [1 1 1 1;  % Reguła 1: wejście MF1 -> wyjście Out1
            2 2 1 1;
            3 3 1 1];  % Reguła 2: wejście MF2 -> wyjście Out2  % Reguła 3: wejście MF3 -> wyjście Out3

% Dodanie reguł do systemu
fis_2 = addRule(fis_2, ruleList);

% Symulacja modelu Hammersteina
U_fuzzy = zeros(obiekt.kk, 1); % Przepuszczenie przez model TS
Y_out = zeros(1, obiekt.kk);

% U = [repelem((rand(1, obiekt.kk/250) * 40 - 20), 250)];
% [Y_real, Y_lin] = obiekt.rk4([U; zeros(1, obiekt.kk)], obiekt.kk); % Symulacja rzeczywistego układu

for k = 1:obiekt.kk
    [~, degrees] = evalfis(fis_2, U(1, k));
    output = 0;
    for i = 1:length(a_param)
        output = output + degrees(i) * (a_param(i)*sinh(U(1, k)/b_param(i)));
    end
    U_fuzzy(k) = output / sum(degrees);

    if k<obiekt.delay+3
        Y_out(k) = 0;
    else
        Y_out(k) = - a*[Y_out(k-1:-1:k-2)]' + b*[U_fuzzy(k-(obiekt.delay+1):-1:k-(obiekt.delay+2))];
    end
end

% Wizualizacja wyników
figure;
plot(t, Y_real, 'b', t, Y_lin, 'g', t, Y_out, 'r');
% legend('RK4', 'RK4 liniowy', 'Hammerstein (optymalny TS)', 'Location', 'southwest');
title('Porównanie wyjścia układu rzeczywistego i modelu');
grid on;

E_lin = sum((Y_real - Y_lin).^2);
E_out = sum((Y_real - Y_out).^2);
fprintf("\nE_lin = %.3f\n", E_lin);
fprintf("E_out = %.3f\n", E_out);

%% Losowość
for j = 1:5
    U = [repelem((rand(1, obiekt.kk/250) * 90 - 45), 250)];
    [Y_real, Y_lin] = obiekt.rk4([U; zeros(1, obiekt.kk)], obiekt.kk); % Symulacja rzeczywistego układu
    
    % Symulacja modelu Hammersteina
    U_fuzzy = zeros(obiekt.kk, 1); % Przepuszczenie przez model TS
    U_fuzzy_2 = zeros(obiekt.kk, 1); % Przepuszczenie przez model TS
    Y_out = zeros(1, obiekt.kk);
    Y_out_2 = zeros(1, obiekt.kk);
    
    for k = obiekt.delay+3:obiekt.kk
        [~, degrees] = evalfis(fis, U(k));
        [~, degrees_2] = evalfis(fis_2, U(k));
        output = 0;
        output_2 = 0;
        for i = 1:length(a_param)
            output = output + degrees(i) * (fis.Outputs.MembershipFunctions(i).Parameters(1)*sinh(U(k)/fis.Outputs.MembershipFunctions(i).Parameters(2)));
            output_2 = output_2 + degrees_2(i) * (fis_2.Outputs.MembershipFunctions(i).Parameters(1)*sinh(U(k)/fis_2.Outputs.MembershipFunctions(i).Parameters(2)));
        end
        U_fuzzy(k) = output / sum(degrees);
        U_fuzzy_2(k) = output_2 / sum(degrees_2);
        Y_out(k) = - a*[Y_out(k-1:-1:k-2)]' + b*[U_fuzzy(k-(obiekt.delay+1):-1:k-(obiekt.delay+2))];
        Y_out_2(k) = - a*[Y_out_2(k-1:-1:k-2)]' + b*[U_fuzzy_2(k-(obiekt.delay+1):-1:k-(obiekt.delay+2))];
    end
    
    % figure;
    % plot(t, Y_real, 'b', t, Y_lin, 'g', t, Y_out, 'r');
    % % legend('RK4', 'RK4 liniowy', 'Hammerstein (optymalny TS)');
    % title('Porównanie wyjścia układu rzeczywistego i modelu');
    % grid on;
    
    E_lin = sum((Y_real - Y_lin).^2);
    E_out = sum((Y_real - Y_out).^2);
    E_out_2 = sum((Y_real - Y_out_2).^2);
    fprintf("\n%d. E_lin = %.3f\n", j, E_lin);
    fprintf("%d. E_out = %.3f\n", j, E_out);
    fprintf("%d. E_out_2 = %.3f\n", j, E_out_2);
end

%% Rysowanie wykresów
u = [ones(1, obiekt.kk)*-15
    ones(1, obiekt.kk)*0
    ones(1, obiekt.kk)*-15];
[y, y_L, E_h, E_pH] = obiekt.modifiedEuler(u, obiekt.kk);

if ~exist('ph_figure', 'var') || ~isvalid(ph_figure)
    ph_figure = figure;
else
    figure(ph_figure); % Jeśli istnieje, przełącz na nią
end

plot(0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp, y(1,:), 'g-', 'LineWidth', 2);
hold on;
plot(0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp, y_L(1,:), 'm-', 'LineWidth', 2);
xlabel('t [s]');
ylabel('h [cm]');
title('Wartość sygnału wyjściowego h - wymuszenie Q_1');
legend('h (nieliniowe)', 'h (liniowe)', 'Location', 'northeast');
grid on;
% saveas(gcf, 'D:/EiTI/MGR/raporty/raport_MGR/pictures/wymuszenie_Q1h.png');  % Zapisuje jako plik PNG

% file = fopen('errors.txt', 'a'); % Otwórz plik do zapisu (tryb 'w' nadpisuje plik)
% fprintf(file, '%.2f\t%.3f\t%.3f\n', u(1,1), E_h, E_pH);
% fclose(file); % Zamknij plik

if ~exist('h_figure', 'var') || ~isvalid(h_figure)
    h_figure = figure;
else
    figure(h_figure); % Jeśli istnieje, przełącz na nią
end

plot(0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp, y(2,:), 'b-', 'LineWidth', 2);
hold on;
plot(0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp, y_L(2,:), 'r-', 'LineWidth', 2);
xlabel('t [s]');
ylabel('pH');
title('Wartość sygnału wyjściowego pH - wymuszenie Q_1');
legend('pH (nieliniowe)', 'pH (liniowe)', 'Location', 'northeast');
grid on;

% saveas(gcf, 'D:/EiTI/MGR/raporty/raport_MGR/pictures/wymuszenie_Q1pH.png');  % Zapisuje jako plik PNG