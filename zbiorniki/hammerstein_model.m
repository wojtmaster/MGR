%% LINEAR HAMMERSTEIN MODEL
clear all;
obiekt = Obiekt();
[~, ~, a, b] = obiekt.linearization(90, 30);

% Zakres sterowania
U_min = -45;
U_max = 45;
U_center = linspace(U_min, U_max, 5); % Środek zbiorów

%Generacja danych sterujących i RK4
U = repelem([0, -22.5, -45, 22.5, 45], 400);
[Y_real, Y_lin] = obiekt.rk4([U; zeros(1, obiekt.kk)], obiekt.kk); % Symulacja rzeczywistego układu
t = (0:length(U)-1) * obiekt.Tp; % Czas w sekundach (próbkowanie = 20s)

%% 2. Tworzenie początkowego systemu rozmytego TS
close all;
fis = sugfis('Name', 'F1_Hammerstein', 'Type', 'sugeno');
fis = addInput(fis, [U_min U_max], 'Name', 'U');

% Definiowanie funkcji przynależności (gaussmf)
for i = 1:length(U_center) 
    fis = addMF(fis, 'U', 'gaussmf', [12, U_center(i)]);
end

% Definiowanie wyjścia i początkowych następników (a_i * u + b_i)
fis = addOutput(fis, [U_min U_max], 'Name', 'F1');

% Początkowe współczynniki (a_i, b_i)
a_param = [0.795 0.9 1.0933 1.07 1.2034];
b_param = [0.0001 0.0002 0.0001 0 0.001];

% Dodanie reguł TS w postaci liniowej
for i = 1:length(U_center)
    fis = addMF(fis, 'F1', 'linear', [a_param(i), b_param(i)]);
end

% Reguły Takagi-Sugeno: [inputMF, outputMF, weight]
ruleList = [1 1 1 1;  % Reguła 1: wejście MF1 -> wyjście Out1
            2 2 1 1;  % Reguła 2: wejście MF2 -> wyjście Out2
            3 3 1 1;
            4 4 1 1;
            5 5 1 1];  % Reguła 3: wejście MF3 -> wyjście Out3

% Dodanie reguł do systemu
fis = addRule(fis, ruleList);

% Symulacja modelu Hammersteina
U_fuzzy = evalfis(fis, U'); % Przepuszczenie przez model TS
Y_out = zeros(1, obiekt.kk);

for k = obiekt.delay+3:obiekt.kk
    Y_out(k) = - a*[Y_out(k-1:-1:k-2)]' + b*[U_fuzzy(k-(obiekt.delay+1):-1:k-(obiekt.delay+2))];
end

% Wizualizacja wyników
figure;
plot(t, Y_real, 'b', t, Y_lin, 'g', t, Y_out, 'r');
legend('RK4', 'RK4 liniowy', 'Hammerstein (optymalny TS)', 'Location', 'southwest');
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
    U_fuzzy = evalfis(fis, U');
    Y_out = zeros(1, obiekt.kk);
    
    for k = obiekt.delay+3:obiekt.kk
        Y_out(k) = - a*[Y_out(k-1:-1:k-2)]' + b*[U_fuzzy(k-(obiekt.delay+1):-1:k-(obiekt.delay+2))];
    end
    
    figure;
    plot(t, Y_real, 'b', t, Y_lin, 'g', t, Y_out, 'r');
    % legend('RK4', 'RK4 liniowy', 'Hammerstein (optymalny TS)');
    title('Porównanie wyjścia układu rzeczywistego i modelu');
    grid on;
    
    E_lin = sum((Y_real - Y_lin).^2);
    E_out = sum((Y_real - Y_out).^2);
    fprintf("\n%d. E_lin = %.3f\n", j, E_lin);
    fprintf("%d. E_out = %.3f\n", j, E_out);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% NONLINEAR HAMMERSTEIN MODEL
clear all;
obiekt = Obiekt();
[~, ~, a, b] = obiekt.linearization(90, 30);

% Zakres sterowania
U_min = -45;
U_max = 45;
U_center = linspace(U_min, U_max, 3); % Środek zbiorów
% U_center = [-30 0 15]; % Środek zbiorów

%Generacja danych sterujących i RK4
U = repelem([0, -12.5, -22.5, -45, 0, 12.5, 22.5, 45], 250);
[Y_real, Y_lin] = obiekt.rk4([U; zeros(1, obiekt.kk)], obiekt.kk); % Symulacja rzeczywistego układu
t = (0:length(U)-1) * obiekt.Tp; % Czas w sekundach (próbkowanie = 20s)


%% FMINSEARCH
fis = sugfis('Name', 'F1_Hammerstein', 'Type', 'sugeno');
fis = addInput(fis, [U_min U_max], 'Name', 'U');

% Definiowanie funkcji przynależności (gaussmf)
fis = addMF(fis, 'U', 'gaussmf', [20, U_center(1)]);
fis = addMF(fis, 'U', 'gaussmf', [20, U_center(2)]);
fis = addMF(fis, 'U', 'gaussmf', [20, U_center(3)]);

% Definiowanie wyjścia i początkowych następników (a_i * u + b_i)
fis = addOutput(fis, [U_min U_max], 'Name', 'F1');

% Początkowe współczynniki (a_i, b_i)
a_param = [76.4472   72.8656  100.0000]; % Można dobrać inaczej
b_param = [100.0000   71.1456   88.1261];

for i = 1:length(U_center)
    fis = addMF(fis, 'F1', 'linear', [a_param(i), b_param(i)]);
end

% Reguły Takagi-Sugeno: [inputMF, outputMF, weight]
ruleList = [1 1 1 1;  % Reguła 1: wejście MF1 -> wyjście Out1
            2 2 1 1;  % Reguła 2: wejście MF2 -> wyjście Out2
            3 3 1 1];  % Reguła 3: wejście MF3 -> wyjście Out3
fis = addRule(fis, ruleList);

% Optymalizacja przy pomocy fminsearch
initial_params = [a_param, b_param];  % Początkowe wartości a_param i b_param
options = optimset('Display', 'iter', 'MaxFunEvals', 1000, 'MaxIter', 1000); % Opcje optymalizacji
optimal_params = fminsearch(@(params) objective_function(params, fis, U, Y_real, obiekt, U_center, a, b), initial_params, options);

% Po optymalizacji
a_optimal = optimal_params(1:length(U_center));
b_optimal = optimal_params(length(U_center)+1:end);

% Wyświetlanie wyników optymalizacji
fprintf('Optymalne parametry a: \n');
disp(a_optimal);
fprintf('Optymalne parametry b: \n');
disp(b_optimal);

%% 2. Tworzenie początkowego systemu rozmytego TS
close all;
fis = sugfis('Name', 'F1_Hammerstein', 'Type', 'sugeno');
fis = addInput(fis, [U_min U_max], 'Name', 'U');

% Definiowanie funkcji przynależności (gaussmf)
fis = addMF(fis, 'U', 'gaussmf', [20, U_center(1)]);
fis = addMF(fis, 'U', 'gaussmf', [20, U_center(2)]);
fis = addMF(fis, 'U', 'gaussmf', [20, U_center(3)]);

% Definiowanie wyjścia i początkowych następników (a_i * u + b_i)
fis = addOutput(fis, [U_min U_max], 'Name', 'F1');

% Początkowe współczynniki (a_i, b_i)
a_param = [60.5 50.7 83.5]; % Można dobrać inaczej
b_param = [80 50 75];

% Dodanie reguł TS w postaci liniowej
for i = 1:length(U_center)
    fis = addMF(fis, 'F1', 'linear', [a_param(i), b_param(i)]);
end

% Reguły Takagi-Sugeno: [inputMF, outputMF, weight]
ruleList = [1 1 1 1;  % Reguła 1: wejście MF1 -> wyjście Out1
            2 2 1 1;
            3 3 1 1];  % Reguła 2: wejście MF2 -> wyjście Out2  % Reguła 3: wejście MF3 -> wyjście Out3

% Dodanie reguł do systemu
fis = addRule(fis, ruleList);

% Symulacja modelu Hammersteina
U_fuzzy = zeros(obiekt.kk, 1); % Przepuszczenie przez model TS
Y_out = zeros(1, obiekt.kk);

% U = [repelem((rand(1, obiekt.kk/250) * 40 - 20), 250)];
% [Y_real, Y_lin] = obiekt.rk4([U; zeros(1, obiekt.kk)], obiekt.kk); % Symulacja rzeczywistego układu

for k = 1:obiekt.kk
    [~, degrees] = evalfis(fis, U(1, k));
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

%% Funkcja celu: oblicza błąd średniokwadratowy E_out
function E_out = objective_function(params, fis, U, Y_real, obiekt, U_center, a, b)
    % Parametry do optymalizacji
    a_param = params(1:length(U_center));  % a_param
    b_param = params(length(U_center)+1:end);  % b_param
    
    % Aktualizacja modelu rozmytego z nowymi parametrami
    for i = 1:length(U_center)
        fis.Outputs.MembershipFunctions(i).Parameters(1) = a_param(i);
        fis.Outputs.MembershipFunctions(i).Parameters(2) = b_param(i);
    end
    
    % Symulacja modelu Hammersteina
    U_fuzzy = zeros(obiekt.kk, 1); % Przepuszczenie przez model TS
    Y_out = zeros(1, obiekt.kk);
    
    for k = 1:obiekt.kk
        [~, degrees] = evalfis(fis, U(1, k));
        output = 0;
        for i = 1:length(a_param)
            output = output + degrees(i) * (a_param(i)*sinh(U(1, k)/b_param(i)));
        end
        U_fuzzy(k) = output / sum(degrees);

        if k < obiekt.delay + 3
            Y_out(k) = 0;
        else
            Y_out(k) = - a * [Y_out(k-1:-1:k-2)]' + b * [U_fuzzy(k - (obiekt.delay+1):-1:k - (obiekt.delay+2))];
        end
    end
    
    % Obliczanie błędu średniokwadratowego
    E_out = sum((Y_real - Y_out).^2); % Błąd MSE
end