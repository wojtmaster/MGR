%% LINEAR HAMMERSTEIN MODEL
clear all;
obiekt = Obiekt();
[a_1, a_2, a_3, b_1, b_2, b_3, G_z] = obiekt.linearization(16.6, 0.55, 15.6);

% Zakres sterowania
U_min = -15;
U_max = 15;
U_center = linspace(U_min, U_max, 5); % Środek zbiorów

for i = 1:30
    if (i <= 15)
        U_tmp = [repelem([0, -7.5, -15, 7.5, 15], 420);
                 zeros(1, obiekt.kk);
                 repelem((rand(1, obiekt.kk/300) * 30 - 15), 300)];
    else
        U_tmp = [repelem((rand(1, obiekt.kk/300) * 30 - 15), 300);
                 zeros(1, obiekt.kk);
                 repelem([0, -7.5, -15, 7.5, 15], 420)];
    end
    [Y_tmp, ~] = obiekt.modifiedEuler(U_tmp, obiekt.kk);

    U_Q1(i, :) = U_tmp(1, :);
    U_Q3(i, :) = U_tmp(3, :);
    Y_train(i, :) = Y_tmp(1, :);
end
t = (0:length(U_Q1(1,:))-1) * obiekt.Tp;

%% 2. Tworzenie początkowego systemu rozmytego TS
close all;
fis = sugfis('Name', 'Hammerstein', 'Type', 'sugeno');
rules_number = 15;

fis = addInput(fis, [U_min U_max], 'Name', 'u1');
fis = addInput(fis, [U_min U_max], 'Name', 'u2');

% Definiowanie funkcji przynależności (gaussmf)
for i = 1:length(U_center) 
    fis = addMF(fis, 'u1', 'gaussmf', [5, U_center(i)]);
    fis = addMF(fis, 'u2', 'gaussmf', [5, U_center(i)]);
end

% Definiowanie wyjścia i początkowych następników (a_i * u + b_i)
fis = addOutput(fis, [U_min U_max], 'Name', 'u_fuzzy');

% Początkowe współczynniki (a_i, b_i, c_i)
a_param = ones(1,rules_number)*0.5;
b_param = ones(1,rules_number)*0.5;
c_param = ones(1,rules_number)*0.5;

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
            1 4 4 1 1;
            4 1 4 1 1;
            1 5 5 1 1;
            5 1 5 1 1;
            2 2 6 1 1;
            2 3 7 1 1;
            3 2 7 1 1;
            2 4 8 1 1;
            4 2 8 1 1;
            2 5 9 1 1;
            5 2 9 1 1;
            3 3 10 1 1;
            3 4 11 1 1;
            4 3 11 1 1;
            3 5 12 1 1;
            5 3 12 1 1;
            4 4 13 1 1;
            4 5 14 1 1;
            5 4 14 1 1;
            5 5 15 1 1];  

% Dodanie reguł do systemu
fis = addRule(fis, ruleList);

% Optymalizacja przy pomocy fminsearch
initial_params = [a_param, b_param c_param];  % Początkowe wartości a_param i b_param
options = optimset('Display', 'iter', 'MaxFunEvals', 1000, 'MaxIter', 1000); % Opcje optymalizacji
optimal_params = fminsearch(@(params) linearCoeff(params, fis, U_Q1, U_Q3, Y_train, rules_number, a_1, b_1.Q_1, b_1.Q_3), initial_params, options);

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

%% Symulacja modelu Hammersteina
for i = 1:rules_number
    fis.Outputs.MembershipFunctions(i).Parameters(1) = a_optimal(i);
    fis.Outputs.MembershipFunctions(i).Parameters(2) = b_optimal(i);
    fis.Outputs.MembershipFunctions(i).Parameters(3) = c_optimal(i);
end

U = [repelem((rand(1, obiekt.kk/300) * 30 - 15), 300);
    zeros(1,obiekt.kk);
    repelem((rand(1, obiekt.kk/300) * 30 - 15), 300)];
[Y_real, Y_lin] = obiekt.modifiedEuler(U, obiekt.kk); % Symulacja rzeczywistego układu

U_fuzzy = evalfis(fis, [U(1,:)', U(3,:)']);
Y_out = zeros(1, obiekt.kk);

for k = 2:obiekt.kk
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

%% Losowość
for j = 1:5
    U = [repelem((rand(1, obiekt.kk/400) * 30 - 15), 400);
        zeros(1,obiekt.kk);
        repelem((rand(1, obiekt.kk/400) * 30 - 15), 400)];
    [Y_real, Y_lin] = obiekt.modifiedEuler(U, obiekt.kk); % Symulacja rzeczywistego układu
    
    % Symulacja modelu Hammersteina
    U_fuzzy = evalfis(fis, [U(1,:)', U(3,:)']);
    Y_out = zeros(1, obiekt.kk);
    
    for k = 2:obiekt.kk
        Y_out(k) = - a_1*Y_out(k-1) + b_1.Q_1 * U_fuzzy(k-1) + b_1.Q_2 * U(2, k-1) + b_1.Q_3 * U_fuzzy(k-1);
    end
    
    figure;
    plot(t, Y_real(1,:), 'b', t, Y_lin(1,:), 'g', t, Y_out, 'r');
    legend('Eulera', 'Euler liniowy', 'Hammerstein (optymalny TS)');
    title('Porównanie wyjścia układu rzeczywistego i modelu');
    grid on;
    
    E_lin = sum((Y_real(1,:) - Y_lin(1,:)).^2);
    E_out = sum((Y_real(1,:) - Y_out).^2);
    fprintf("\n%d. E_lin = %.3f\n", j, E_lin);
    fprintf("%d. E_out = %.3f\n", j, E_out);
end

%% Funkcja do optymalizacji współczynników
function E_out = linearCoeff(params, fis, U_Q1, U_Q3, Y_real, rules_number, a1, b1_q1, b1_q3)
    % Parametry do optymalizacji
    a_param = params(1:rules_number);
    b_param = params(rules_number+1:2*rules_number);
    c_param = params(2*rules_number+1:end);

    % Zaktualizuj model fuzzy z nowymi parametrami
    for i = 1:rules_number
        fis.Outputs.MembershipFunctions(i).Parameters(1) = a_param(i);
        fis.Outputs.MembershipFunctions(i).Parameters(2) = b_param(i);
        fis.Outputs.MembershipFunctions(i).Parameters(3) = c_param(i);
    end

    % Inicjalizacja błędu
    [M, N] = size(U_Q1);  % M = liczba trajektorii
    E_out = 0;

    for traj = 1:M
        u_q1 = U_Q1(traj,:);
        u_q3 = U_Q3(traj,:);
        y_real_single = Y_real(traj, :);  % Odpowiadające wyjście

        u_fuzzy = evalfis(fis, [u_q1', u_q3']);  % N × 2
        y_out = zeros(1, N);

        for k = 2:N
            y_out(k) = -a1 * y_out(k-1) + b1_q1 * u_fuzzy(k-1) + b1_q3 * u_fuzzy(k-1);
        end

        % Sumuj błąd dla tej trajektorii
        E_out = E_out + sum((y_real_single - y_out).^2);
    end
end