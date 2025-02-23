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
a_param = [77 69 65.2]; % Można dobrać inaczej
b_param = [100 67 60];

% % Początkowe współczynniki (a_i, b_i)
% a_param = [74.5 78 50.5]; % Można dobrać inaczej
% b_param = [100 70 49];

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

for k = obiekt.delay+3:obiekt.kk
    [~, degrees] = evalfis(fis, U(k));
    output = 0;
    for i = 1:length(a_param)
        output = output + degrees(i) * (a_param(i)*sinh(U(k)/b_param(i)));
    end
    U_fuzzy(k) = output / sum(degrees);
    Y_out(k) = - a*[Y_out(k-1:-1:k-2)]' + b*[U_fuzzy(k-(obiekt.delay+1):-1:k-(obiekt.delay+2))];
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
    Y_out = zeros(1, obiekt.kk);
    
    for k = obiekt.delay+3:obiekt.kk
        [~, degrees] = evalfis(fis, U(k));
        output = 0;
        for i = 1:length(a_param)
            output = output + degrees(i) * (a_param(i)*sinh(U(k)/b_param(i)));
        end
        U_fuzzy(k) = output / sum(degrees);
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