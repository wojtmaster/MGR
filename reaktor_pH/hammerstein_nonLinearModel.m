%% NONLINEAR HAMMERSTEIN MODEL
clear all;
obiekt = Obiekt();
[a_1, a_2, a_3, b_1, b_2, b_3, G_z] = obiekt.linearization(16.6, 0.55, 15.6);

% Zakres sterowania
U_min = -15;
U_max = 15;
U_center = linspace(U_min, U_max, 3); % Środek zbiorów
sigma = 8;

for i = 1:30
    if (i <= 15)
        U_tmp = [repelem([0, -7.5, -15, 7.5, 15], 1000);
                 zeros(1, obiekt.kk);
                 repelem((rand(1, obiekt.kk/500) * (i*2) - i), 500)];
    else
        U_tmp = [repelem((rand(1, obiekt.kk/500) * (i*2-30) - (i-15)), 500);
                zeros(1, obiekt.kk);
                repelem([0, -7.5, -15, 7.5, 15], 1000)];
    end
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
    fis = addMF(fis, 'u1', 'gaussmf', [sigma, U_center(i)]);
    fis = addMF(fis, 'u2', 'gaussmf', [sigma, U_center(i)]);
end

% Definiowanie wyjścia i początkowych następników (a_i * u + b_i)
fis = addOutput(fis, [U_min U_max], 'Name', 'u_fuzzy');

% Początkowe współczynniki (a_i, b_i, c_i)
a_param = ones(1,rules_number)*19;
b_param = ones(1,rules_number)*40;
c_param = ones(1,rules_number)*19;
d_param = ones(1,rules_number)*40;
% a_param = rand(1,rules_number)*1;
% b_param = rand(1,rules_number)*25;
% c_param = rand(1,rules_number)*1;
% d_param = rand(1,rules_number)*25;
% a_param = [1.0711    1.3295   -0.1572    1.5354    1.1508    0.7317];
% b_param = [5.8377    6.2246    5.5145    5.1008    5.0795   -1.0823];
% c_param = [1.7000    1.2547   -0.1612    1.1328    1.9847    0.5973];
% d_param = [9.4321    6.0064   10.7646    4.7766    6.1933    0.3843];

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
initial_params = [a_param, b_param c_param d_param];  % Początkowe wartości a_param i b_param
options = optimset('Display', 'iter', 'MaxFunEvals', 1500, 'MaxIter', 1000); % Opcje optymalizacji
optimal_params = fminsearch(@(params) nonlinearCoeff(params, fis, U_Q1, U_Q3, Y_train, U_center, sigma, rules_number, a_1, b_1.Q_1, b_1.Q_3), initial_params, options);

% Po optymalizacji
a_optimal = optimal_params(1:rules_number);
b_optimal = optimal_params(rules_number+1:2*rules_number);
c_optimal = optimal_params(2*rules_number+1:3*rules_number);
d_optimal = optimal_params(3*rules_number+1:end);

% Wyświetlanie wyników optymalizacji
fprintf('Optymalne parametry a: \n');
disp(a_optimal);
fprintf('Optymalne parametry b: \n');
disp(b_optimal);
fprintf('Optymalne parametry c: \n');
disp(c_optimal);
fprintf('Optymalne parametry d: \n');
disp(d_optimal);

%% 2. Tworzenie początkowego systemu rozmytego TS
gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));

for i = 1:rules_number
    fis.Outputs.MembershipFunctions(i).Parameters(1) = a_optimal(i);
    fis.Outputs.MembershipFunctions(i).Parameters(2) = b_optimal(i);
    fis.Outputs.MembershipFunctions(i).Parameters(3) = c_optimal(i);
end

% Symulacja modelu Hammersteina
U_fuzzy = zeros(1, obiekt.kk); % Przepuszczenie przez model TS
Y_out = zeros(1, obiekt.kk);
rule_to_param = [1 2 2 3 3 1 5 5 4];

U = [repelem((rand(1, obiekt.kk/500) * 30 - 15), 500);
    zeros(1, obiekt.kk);
    repelem((rand(1, obiekt.kk/500) * 30 - 15), 500)];
[Y_real, Y_lin] = obiekt.modifiedEuler(U, obiekt.kk); % Symulacja rzeczywistego układu

for k = 1:obiekt.kk         
    deg_u1 = [gaussmf_val(U(1,k), sigma, U_center(1)), gaussmf_val(U(1,k), sigma, U_center(2)), gaussmf_val(U(1,k), sigma, U_center(3))];
    deg_u2 = [gaussmf_val(U(3,k), sigma, U_center(1)), gaussmf_val(U(3,k), sigma, U_center(2)), gaussmf_val(U(3,k), sigma, U_center(3))];
    
    % Teraz tworzysz reguły ręcznie:
    % Np. reguła 1: u1_1 & u2_1
    degrees_all(k, 1) = deg_u1(1) * deg_u2(1);
    degrees_all(k, 2) = deg_u1(1) * deg_u2(2);
    degrees_all(k, 3) = deg_u1(2) * deg_u2(1);
    degrees_all(k, 4) = deg_u1(1) * deg_u2(3);
    degrees_all(k, 5) = deg_u1(3) * deg_u2(1);
    degrees_all(k, 6) = deg_u1(2) * deg_u2(2);
    degrees_all(k, 7) = deg_u1(2) * deg_u2(3);
    degrees_all(k, 8) = deg_u1(3) * deg_u2(2);
    degrees_all(k, 9) = deg_u1(3) * deg_u2(3);

    % degrees_all(k, 1) = min(deg_u1(1), deg_u2(1));
    % degrees_all(k, 2) = min(deg_u1(1), deg_u2(2));
    % degrees_all(k, 3) = min(deg_u1(2), deg_u2(1));
    % degrees_all(k, 4) = min(deg_u1(1), deg_u2(3));
    % degrees_all(k, 5) = min(deg_u1(3), deg_u2(1));
    % degrees_all(k, 6) = min(deg_u1(2), deg_u2(2));
    % degrees_all(k, 7) = min(deg_u1(2), deg_u2(3));
    % degrees_all(k, 8) = min(deg_u1(3), deg_u2(2));
    % degrees_all(k, 9) = min(deg_u1(3), deg_u2(3));
end

for k = 2:obiekt.kk
    [~, degrees] = evalfis(fis, [U(1,k), U(3,k)]);
    output = 0;
    w = 0;
    for i = 1:length(fis.Rules)
        index = rule_to_param(i);
        % output = output + degrees_all(k,i)*(a_optimal(index)*sinh(U(1,k)/b_optimal(index)) + ...
        %     c_optimal(index)*sinh(U(3,k)/d_optimal(index)));
        
        output = output + degrees_all(k,i)*(a_optimal(index)*sinh(U(1,k)/b_optimal(index)) + ...
            c_optimal(index)*sinh(U(3,k)/d_optimal(index)));
        % output = output + degrees_all(k,i)*(a_optimal(index)*sinh((U(1,k)+U(3,k))/b_optimal(index)));
        w = w + degrees_all(k,i);
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

E_lin = sum((Y_real(1,:) - Y_lin(1,:)).^2) / obiekt.kk;
E_out = sum((Y_real(1,:) - Y_out).^2) / obiekt.kk;
fprintf("\nE_lin = %.3f\n", E_lin);
fprintf("E_out = %.3f\n", E_out);

%% Losowość
close all;
gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));
rule_to_param = [1 2 2 3 3 1 5 5 4];
for j = 1:5
    % Symulacja modelu Hammersteina
    U_fuzzy = zeros(1, obiekt.kk); % Przepuszczenie przez model TS
    Y_out = zeros(1, obiekt.kk);
    
    U = [repelem((rand(1, obiekt.kk/500) * 30 - 15), 500);
        zeros(1, obiekt.kk);
        repelem((rand(1, obiekt.kk/500) * 30 - 15), 500)];
    [Y_real, Y_lin] = obiekt.modifiedEuler(U, obiekt.kk); % Symulacja rzeczywistego układu
    
    for k = 1:obiekt.kk         
        deg_u1 = [gaussmf_val(U(1,k), sigma, U_center(1)), gaussmf_val(U(1,k), sigma, U_center(2)), gaussmf_val(U(1,k), sigma, U_center(3))];
        deg_u2 = [gaussmf_val(U(3,k), sigma, U_center(1)), gaussmf_val(U(3,k), sigma, U_center(2)), gaussmf_val(U(3,k), sigma, U_center(3))];
        
        % Teraz tworzysz reguły ręcznie:
        % Np. reguła 1: u1_1 & u2_1
        degrees_all(k, 1) = deg_u1(1) * deg_u2(1);
        degrees_all(k, 2) = deg_u1(1) * deg_u2(2);
        degrees_all(k, 3) = deg_u1(2) * deg_u2(1);
        degrees_all(k, 4) = deg_u1(1) * deg_u2(3);
        degrees_all(k, 5) = deg_u1(3) * deg_u2(1);
        degrees_all(k, 6) = deg_u1(2) * deg_u2(2);
        degrees_all(k, 7) = deg_u1(2) * deg_u2(3);
        degrees_all(k, 8) = deg_u1(3) * deg_u2(2);
        degrees_all(k, 9) = deg_u1(3) * deg_u2(3);
    
        % degrees_all(k, 1) = min(deg_u1(1), deg_u2(1));
        % degrees_all(k, 2) = min(deg_u1(1), deg_u2(2));
        % degrees_all(k, 3) = min(deg_u1(2), deg_u2(1));
        % degrees_all(k, 4) = min(deg_u1(1), deg_u2(3));
        % degrees_all(k, 5) = min(deg_u1(3), deg_u2(1));
        % degrees_all(k, 6) = min(deg_u1(2), deg_u2(2));
        % degrees_all(k, 7) = min(deg_u1(2), deg_u2(3));
        % degrees_all(k, 8) = min(deg_u1(3), deg_u2(2));
        % degrees_all(k, 9) = min(deg_u1(3), deg_u2(3));
    end
    
    for k = 2:obiekt.kk
        [~, degrees] = evalfis(fis, [U(1,k), U(3,k)]);
        output = 0;
        w = 0;
        for i = 1:length(fis.Rules)
            index = rule_to_param(i);
            output = output + degrees_all(k,i)*(a_optimal(index)*sinh(U(1,k)/b_optimal(index)) + ...
                c_optimal(index)*sinh(U(3,k)/d_optimal(index)));
            % output = output + degrees_all(k,i)*(a_optimal(index)*sinh((U(1,k) + U(3,k)) / b_optimal(index)));
            w = w + degrees_all(k,i);
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
    
    E_lin = sum((Y_real(1,:) - Y_lin(1,:)).^2) / obiekt.kk;
    E_out = sum((Y_real(1,:) - Y_out).^2) / obiekt.kk;
    fprintf("\n%d. E_lin = %.3f\n", j, E_lin);
    fprintf("%d. E_out = %.3f\n", j, E_out);
end


%% Funkcja optymalizująca parametry
function E_out = nonlinearCoeff(params, fis, U_Q1, U_Q3, Y_real, U_center, sigma, rules_number, a1, b1_q1, b1_q3)

    gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));

    % Parametry do optymalizacji
    a_param = params(1:rules_number);
    b_param = params(rules_number+1:2*rules_number);
    c_param = params(2*rules_number+1:end);
    d_param = params(3*rules_number+1:end);

    rule_to_param = [1 2 2 3 3 1 5 5 4];

    % Zaktualizuj model fuzzy z nowymi parametrami
    for i = 1:rules_number
        fis.Outputs.MembershipFunctions(i).Parameters(1) = a_param(i);
        fis.Outputs.MembershipFunctions(i).Parameters(2) = b_param(i);
        % fis.Outputs.MembershipFunctions(i).Parameters(3) = c_param(i);
    end

    % Inicjalizacja błędu
    [M, N] = size(U_Q1);  % M = liczba trajektorii
    degrees_all = zeros(N, length(fis.Rules));
    E_out = 0;

    for traj = 1:M
        u_q1 = U_Q1(traj,:);
        u_q3 = U_Q3(traj,:);
        y_real_single = Y_real(traj, :);  % Odpowiadające wyjście

        for k = 1:N            
            deg_u1 = [gaussmf_val(u_q1(k), sigma, U_center(1)), gaussmf_val(u_q1(k), sigma, U_center(2)), gaussmf_val(u_q1(k), sigma, U_center(3))];
            deg_u2 = [gaussmf_val(u_q3(k), sigma, U_center(1)), gaussmf_val(u_q3(k), sigma, U_center(2)), gaussmf_val(u_q3(k), sigma, U_center(3))];
         
            % Teraz tworzysz reguły ręcznie:
            % Np. reguła 1: u1_1 & u2_1
            degrees_all(k, 1) = deg_u1(1) * deg_u2(1);
            degrees_all(k, 2) = deg_u1(1) * deg_u2(2);
            degrees_all(k, 3) = deg_u1(2) * deg_u2(1);
            degrees_all(k, 4) = deg_u1(1) * deg_u2(3);
            degrees_all(k, 5) = deg_u1(3) * deg_u2(1);
            degrees_all(k, 6) = deg_u1(2) * deg_u2(2);
            degrees_all(k, 7) = deg_u1(2) * deg_u2(3);
            degrees_all(k, 8) = deg_u1(3) * deg_u2(2);
            degrees_all(k, 9) = deg_u1(3) * deg_u2(3);
        % 
        %     % degrees_all(k, 1) = min(deg_u1(1), deg_u2(1));
        %     % degrees_all(k, 2) = min(deg_u1(1), deg_u2(2));
        %     % degrees_all(k, 3) = min(deg_u1(2), deg_u2(1));
        %     % degrees_all(k, 4) = min(deg_u1(1), deg_u2(3));
        %     % degrees_all(k, 5) = min(deg_u1(3), deg_u2(1));
        %     % degrees_all(k, 6) = min(deg_u1(2), deg_u2(2));
        %     % degrees_all(k, 7) = min(deg_u1(2), deg_u2(3));
        %     % degrees_all(k, 8) = min(deg_u1(3), deg_u2(2));
        %     % degrees_all(k, 9) = min(deg_u1(3), deg_u2(3));
        end

        y_out = zeros(1, N);
        u_fuzzy = zeros(1,N);
        for k = 2:N
            % [~, degrees] = evalfis(fis, [u_q1(k), u_q3(k)]);
            output = 0;
            w = 0;
            for i = 1:length(fis.Rules)
                index = rule_to_param(i);
                output = output +  degrees_all(k,i)* (a_param(index)*sinh(u_q1(k)/b_param(index)) + ...
                    c_param(index)*sinh(u_q3(k)/d_param(index)));

                %%% Najlepsze do tej pory
                % output = output +  degrees_all(k,i) * (a_param(index)*sinh((u_q1(k)+u_q3(k))/b_param(index)));
                w = w + degrees_all(k,i);


                % w_curr = degrees(i,1)*degrees(i,2);
                % output = output +  w_curr * (a_param(index)*sinh(u_q1(k)/b_param(index)) + ...
                %     c_param(index)*sinh(u_q3(k)/d_param(index)));
                % w = w + w_curr;
            end
            u_fuzzy(k) = output / w;
            y_out(k) = - a1*y_out(k-1) + b1_q1*u_fuzzy(k-1) + b1_q3*u_fuzzy(k-1);
        end

        % Sumuj błąd dla tej trajektorii
        E_out = E_out + sum((y_real_single - y_out).^2);
    end
end