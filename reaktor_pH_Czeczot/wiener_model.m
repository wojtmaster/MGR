%% LINEAR WIENER MODEL

clear all;
clc;
obiekt = Obiekt();
[~, ~, a, b] = obiekt.linearization(90, 30);

% Zakres wyjścia liniowego bloku dynamicznego
Y_min = -27;
Y_max = 27;
Y_center = linspace(Y_min, Y_max, 5); % Środki zbiorów rozmytych

% Generacja danych sterujących i RK4
U = repelem([0, -22.5, -45, 22.5, 45], 400);
[Y_real, Y_lin] = obiekt.rk4([U; zeros(1, obiekt.kk)], obiekt.kk); % Symulacja rzeczywistego układu
t = (0:length(U)-1) * obiekt.Tp; % Czas w sekundach (próbkowanie = 20s)

%% Tworzenie systemu rozmytego TS dla Wienera
fis = sugfis('Name', 'F1_Wiener', 'Type', 'sugeno');
fis = addInput(fis, [Y_min Y_max], 'Name', 'Y_lin');

% Definiowanie funkcji przynależności (gaussmf)
for i = 1:length(Y_center) 
    fis = addMF(fis, 'Y_lin', 'gaussmf', [10, Y_center(i)]);
end

% Definiowanie wyjścia i początkowych następników (a_i * y + b_i)
fis = addOutput(fis, [Y_min Y_max], 'Name', 'F1');

% Początkowe współczynniki (a_i, b_i)
a_param = [0.77 0.92 1 1.05 1.24]; % Można dostroić
b_param = [0 0 0 0 0];

% Dodanie reguł TS w postaci liniowej
for i = 1:length(Y_center)
    fis = addMF(fis, 'F1', 'linear', [a_param(i), b_param(i)]);
end

% Reguły Takagi-Sugeno: [inputMF, outputMF, weight]
ruleList = [1 1 1 1;
            2 2 1 1;
            3 3 1 1;
            4 4 1 1;
            5 5 1 1];

% Dodanie reguł do systemu
fis = addRule(fis, ruleList);

% Symulacja modelu Wienera
y_out = zeros(1, obiekt.kk);
Y_out = zeros(1, obiekt.kk);
for k = obiekt.delay+3:obiekt.kk
    y_out(k) = - a*[y_out(k-1:-1:k-2)]' + b*[U(k-(obiekt.delay+1):-1:k-(obiekt.delay+2))]';
    Y_out(k) = evalfis(fis, y_out(k));
end

% Wizualizacja wyników
figure;
plot(t, Y_real, 'b', t, Y_lin, 'g', t, Y_out, 'r');
% legend('RK4', 'RK4 liniowy', 'Wiener (optymalny TS)', 'Location', 'northwest');
title('Porównanie wyjścia układu rzeczywistego i modelu');

E_lin = sum((Y_real - Y_lin).^2);
E_out = sum((Y_real - Y_out).^2);
fprintf("\nE_lin = %.3f\n", E_lin);
fprintf("E_out = %.3f\n", E_out);

%%
for j = 1:5
    U = [repelem((rand(1, obiekt.kk/250) * 90 - 45), 250)];
    [Y_real, Y_lin] = obiekt.rk4([U; zeros(1, obiekt.kk)], obiekt.kk); % Symulacja rzeczywistego układu
    
    % Symulacja modelu Wienera
    y_out = zeros(1, obiekt.kk);
    Y_out = zeros(1, obiekt.kk);
    for k = obiekt.delay+3:obiekt.kk
        y_out(k) = - a*[y_out(k-1:-1:k-2)]' + b*[U(k-(obiekt.delay+1):-1:k-(obiekt.delay+2))]';
        Y_out(k) = evalfis(fis, y_out(k));
    end
    
    figure;
    plot(t, Y_real, 'b', t, Y_lin, 'g', t, Y_out, 'r');
    % legend('RK4', 'RK4 liniowy', 'Wiener (optymalny TS)', 'Location', 'northwest');
    title('Porównanie wyjścia układu rzeczywistego i modelu');
    
    E_lin = sum((Y_real - Y_lin).^2);
    E_out = sum((Y_real - Y_out).^2);
    fprintf("\n%d. E_lin = %.3f\n", j, E_lin);
    fprintf("%d. E_out = %.3f\n", j, E_out);
    grid on;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% NONLINEAR WIENER MODEL
clear all;
clc;
obiekt = Obiekt();
[~, ~, a, b] = obiekt.linearization(90, 30);

% Zakres wyjścia liniowego bloku dynamicznego
Y_min = -27;
Y_max = 27;
Y_center = linspace(Y_min, Y_max, 3); % Środki zbiorów rozmytych

% Generacja danych sterujących i RK4
U = repelem([0, -12.5, -22.5, -45, 0, 12.5, 22.5, 45], 250);
[Y_real, Y_lin] = obiekt.rk4([U; zeros(1, obiekt.kk)], obiekt.kk); % Symulacja rzeczywistego układu
t = (0:length(U)-1) * obiekt.Tp; % Czas w sekundach (próbkowanie = 20s)

%% Tworzenie systemu rozmytego TS dla Wienera
fis = sugfis('Name', 'F1_Wiener', 'Type', 'sugeno');
fis = addInput(fis, [Y_min Y_max], 'Name', 'Y_lin');

% Definiowanie funkcji przynależności (gaussmf)
for i = 1:length(Y_center) 
    fis = addMF(fis, 'Y_lin', 'gaussmf', [15, Y_center(i)]);
end

% Definiowanie wyjścia i początkowych następników (a_i * y + b_i)
fis = addOutput(fis, [Y_min Y_max], 'Name', 'F1');

% Początkowe współczynniki (a_i, b_i)
a_param = [44.1 47.2 50.7]; % Można dostroić
b_param = [60 45 45];

% Dodanie reguł TS w postaci liniowej
for i = 1:length(Y_center)
    fis = addMF(fis, 'F1', 'linear', [a_param(i), b_param(i)]);
end

% Reguły Takagi-Sugeno: [inputMF, outputMF, weight]
ruleList = [1 1 1 1;
            2 2 1 1;
            3 3 1 1];

% Dodanie reguł do systemu
fis = addRule(fis, ruleList);

% Symulacja modelu Wienera
y_out = zeros(1, obiekt.kk);
Y_out = zeros(1, obiekt.kk);
for k = obiekt.delay+3:obiekt.kk
    y_out(k) = - a*[y_out(k-1:-1:k-2)]' + b*[U(1, k-(obiekt.delay+1):-1:k-(obiekt.delay+2))]';
    % y_out(k) = max(Y_min, min(Y_max, y_out(k)));
    [~, degrees] = evalfis(fis, y_out(k)); % Przepuszczenie przez model TS
    output = 0;
    for i = 1:length(a_param)
        output = output + degrees(i) * (a_param(i) * sinh(y_out(k) / b_param(i)));
    end
    Y_out(k) = output / sum(degrees);
    % Y_out(k) = max(Y_min, min(Y_max, Y_out(k)));
end

% Wizualizacja wyników
figure;
plot(t, Y_real, 'b', t, Y_lin, 'g', t, Y_out, 'r');
% legend('RK4', 'RK4 liniowy', 'Wiener (optymalny TS)', 'Location', 'northwest');
title('Porównanie wyjścia układu rzeczywistego i modelu');

E_lin = sum((Y_real - Y_lin).^2);
E_out = sum((Y_real - Y_out).^2);
fprintf("\nE_lin = %.3f\n", E_lin);
fprintf("E_out = %.3f\n", E_out);

%% Losowość
for j = 1:5
    U = [repelem((rand(1, obiekt.kk/250) * 90 - 45), 250)];
    [Y_real, Y_lin] = obiekt.rk4([U; zeros(1, obiekt.kk)], obiekt.kk); % Symulacja rzeczywistego układu
    
    % Symulacja modelu Wienera
    y_out = zeros(1, obiekt.kk);
    Y_out = zeros(1, obiekt.kk);
    for k = obiekt.delay+3:obiekt.kk
        y_out(k) = - a*[y_out(k-1:-1:k-2)]' + b*[U(k-(obiekt.delay+1):-1:k-(obiekt.delay+2))]';
        % y_out(k) = max(Y_min, min(Y_max, y_out(k)));
        [~, degrees] = evalfis(fis, y_out(k)); % Przepuszczenie przez model TS
        output = 0;
        for i = 1:length(a_param)
            output = output + degrees(i) * (a_param(i) * sinh(y_out(k) / b_param(i)));
        end
        Y_out(k) = output / sum(degrees);
        % Y_out(k) = max(Y_min, min(Y_max, Y_out(k)));
    end
    
    figure;
    plot(t, Y_real, 'b', t, Y_lin, 'g', t, Y_out, 'r');
    % legend('RK4', 'RK4 liniowy', 'Wiener (optymalny TS)', 'Location', 'northwest');
    title('Porównanie wyjścia układu rzeczywistego i modelu');
    
    E_lin = sum((Y_real - Y_lin).^2);
    E_out = sum((Y_real - Y_out).^2);
    fprintf("\n%d. E_lin = %.3f\n", j, E_lin);
    fprintf("%d. E_out = %.3f\n", j, E_out);
end