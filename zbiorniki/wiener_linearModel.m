%% LINEAR HAMMERSTEIN MODEL
clear all;
obiekt = Obiekt();
obiekt.linearization(90, 30);

%% Identyfikacja dynamiki
U = [ones(1, 100)
    zeros(1, 100)];
[Y_step, ~] = obiekt.modifiedEuler(U, 100);
% Normalizacja odpowiedzi skokowej
Y_step = Y_step / Y_step(end);
t = (0:length(U)-1) * obiekt.Tp; % Czas w sekundach (próbkowanie = 20s)

% Funkcja celu (error)
fun = @(T) model_error(T, U, Y_step, obiekt);

% Startowa wartość (70)
% T0 = [50, 50];
% Optymalna stała czasowa 
T0 = [193.4125      36.78401];

% Szukanie optymalnej stałej czasowej
% options = optimset('Display', 'iter', 'MaxFunEvals', 2000, 'MaxIter', 1000); % Opcje optymalizacji
% T_opt = fminsearch(fun, T0, options);
T_opt = T0;

% Wyświetlenie wyniku
disp(['Optymalna stała czasowa T = ', num2str(T_opt)]);

G = tf(1, conv([T_opt(1), 1], [T_opt(2), 1]));
G.InputDelay = obiekt.tau;
G_z = c2d(G, obiekt.Tp, 'zoh');
G_z.Variable = 'z^-1';

Y = zeros(size(Y_step));
Y(1:7) = Y_step(1:7);
for k = 8:length(U)
    Y(k) = - G_z.Denominator{1}(2)*Y(k-1) - G_z.Denominator{1}(3)*Y(k-2) ...
        + G_z.Numerator{1}(2)*U(1, k-6) + G_z.Numerator{1}(3)*U(1, k-7) ...
        + G_z.Numerator{1}(2)*U(2, k-1) + G_z.Numerator{1}(3)*U(2, k-2);
end

plot(t, Y_step, 'b', t, Y, 'g');
grid on;

%% Tworzenie początkowego systemu rozmytego TS
close all;

%Generacja danych sterujących i RK4
U = [linspace(-45, 45, 100);
    linspace(-15, 15, 100)];
[F1_grid, FD_grid] = meshgrid(U(1,:), U(2,:));  % kombinacje
Y = ((obiekt.F_10+F1_grid + obiekt.F_D0+FD_grid) / obiekt.alpha_2).^2 - obiekt.h_20;
t = (0:length(U)-1) * obiekt.Tp; % Czas w sekundach (próbkowanie = 20s)

h = zeros(100, 100);
y = zeros(1, 100);
for i = 1:100
    for j = 1:100
        u = [ones(1,100) * U(1,j);
            ones(1,100) * U(2,i)];
        for k = 8:length(u)
            y(k) = - G_z.Denominator{1}(2)*y(k-1) - G_z.Denominator{1}(3)*y(k-2) ...
                + G_z.Numerator{1}(2)*u(1, k-6) + G_z.Numerator{1}(3)*u(1, k-7) ...
                + G_z.Numerator{1}(2)*u(2, k-1) + G_z.Numerator{1}(3)*u(2, k-2);
        end
        % h(i,j) = ((120+y(end))/20)^2 - 36;
        h(i,j) = ((90+u(1, k-6) + 30 + u(2, k-1))/20)^2 - 36;
    end
end

figure;
surf(F1_grid, FD_grid, h);
figure;
surf(F1_grid, FD_grid, Y);

%% System rozmyty
% Zakres sterowania
Y_min = min(min(h));
Y_max = max(max(h));
Y_center = linspace(Y_min, Y_max, 5); % Środek zbiorów
sigma = 15;

fis = sugfis('Name', 'Wiener', 'Type', 'sugeno');
fis = addInput(fis, [Y_min Y_max], 'Name', 'Y');

% Definiowanie funkcji przynależności (gaussmf)
for i = 1:length(Y_center) 
    fis = addMF(fis, 'Y', 'gaussmf', [15, Y_center(i)]);
end

% Definiowanie wyjścia i początkowych następników (a_i * u + b_i)
fis = addOutput(fis, [Y_min Y_max], 'Name', 'Y_fuzzy');

% % Początkowe współczynniki (a_i, b_i)
% a_param = ones(1,5)*0.6;
% b_param = ones(1,5)*1.9;

a_param = [	0.5079	0.5508	0.5985	0.6438	0.6801];
b_param = [	3.5118	0.5148	-0.1804	0.6402	4.2028];

% a_param = [0.5409    0.5918    0.6604    0.7429    0.8309];
% b_param = [41.5737   37.8453   36.0261   34.0338   31.3233];

% Dodanie reguł TS w postaci liniowej
for i = 1:length(Y_center)
    fis = addMF(fis, 'Y_fuzzy', 'linear', [a_param(i), b_param(i)]);
end

% Reguły Takagi-Sugeno: [inputMF, outputMF, weight]
ruleList = [1 1 1 1;
            2 2 1 1;
            3 3 1 1;
            4 4 1 1;
            5 5 1 1];

% Dodanie reguł do systemu
fis = addRule(fis, ruleList);

% % Optymalizacja przy pomocy fminsearch
% initial_params = [a_param, b_param];  % Początkowe wartości a_param i b_param
% options = optimset('Display', 'iter', 'MaxFunEvals', 2000, 'MaxIter', 1000); % Opcje optymalizacji
% optimal_params = fminsearch(@(params) linearCoeff(params, fis, sigma, Y_center, h, Y, length(fis.Rules)), initial_params, options);

% % Po optymalizacji
% a_optimal = optimal_params(1:length(fis.Rules));
% b_optimal = optimal_params(length(fis.Rules)+1:end);
a_optimal = a_param;
b_optimal = b_param;

% Wyświetlanie wyników optymalizacji
fprintf('a_param = [');
fprintf("\t%.4f", a_optimal);
fprintf('];\n');
fprintf('b_param = [');
fprintf("\t%.4f", b_optimal);
fprintf('];\n');

%% Check fuzzy static
for i = 1:length(fis.Rules)
    fis.Outputs.MembershipFunctions(i).Parameters(1) = a_optimal(i);
    fis.Outputs.MembershipFunctions(i).Parameters(2) = b_optimal(i);
end

Y_fuzzy = zeros(size(h));
for i = 1:length(h)
    for j = 1:length(h)
        Y_fuzzy(i,j) = evalfis(fis, h(i,j)) + obiekt.h_20;
    end
end

surf(F1_grid, FD_grid, Y_fuzzy);
% plot(h, Y, 'b-', h, Y_fuzzy, 'ro');
% grid on;

%% Losowość
U = [repelem((rand(1, obiekt.kk/400) * 90 - 45), 400)
    repelem((rand(1, obiekt.kk/500) * 30 - 15), 500)];
[Y_real, Y_lin] = obiekt.rk4(U, obiekt.kk);
t = (0:length(U)-1) * obiekt.Tp; % Czas w sekundach (próbkowanie = 20s)

Y_0 = evalfis(fis, 0);
Y = zeros(1,obiekt.kk);
Y_fuzzy = zeros(1, obiekt.kk);

for k = 8:obiekt.kk
    Y(k) = - G_z.Denominator{1}(2)*Y(k-1) - G_z.Denominator{1}(3)*Y(k-2) ...
        + G_z.Numerator{1}(2)*U(1, k-6) + G_z.Numerator{1}(3)*U(1, k-7) ...
        + G_z.Numerator{1}(2)*U(2, k-1) + G_z.Numerator{1}(3)*U(2, k-2);
    Y_fuzzy(k) = evalfis(fis, Y(k));
end

plot(t, Y_real, 'b');
hold on;
plot(t, Y_lin, 'r');
plot(t, Y_fuzzy, 'g');
grid on;
legend('Y_{real}', 'Y_{lin}', 'Y_fuzzy');

%% Funkcja do optymalizacji współczynników
function E_out = linearCoeff(params, fis, sigma, center, U, Y, rules_number)
    gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));

    % Parametry do optymalizacji
    a_param = params(1:rules_number);
    b_param = params(rules_number+1:end);

    for i = 1:rules_number
        fis.Outputs.MembershipFunctions(i).Parameters(1) = a_param(i);
        fis.Outputs.MembershipFunctions(i).Parameters(2) = b_param(i);
    end
    
    Y_out = zeros(size(U));
    E_out = 0;
    % Oblicz odpowiedzi systemu rozmytego dla wszystkich U
    for i = 1:length(U)
        for j = 1:length(U)
            degrees = [gaussmf_val(U(i,j), sigma, center(1)), gaussmf_val(U(i,j), sigma, center(2)), ...
                gaussmf_val(U(i,j), sigma, center(3)), gaussmf_val(U(i,j), sigma, center(4)), gaussmf_val(U(i,j), sigma, center(5))];
            w = 0;
            output = 0;
            for k = 1:rules_number
                output = output + degrees(k)*(a_param(k)*U(i,j) + b_param(k));
                w = w + degrees(k);
            end
            Y_out(i,j) = output / w;
            % Y_out(i,j) = evalfis(fis, U(i,j));
            E_out = E_out + sum((Y(i,j) - Y_out(i,j))^2);
        end
    end
end

function E = model_error(T, U, Y_step, obiekt)
    % Tworzenie nowej transmitancji z aktualnym T
    G = tf(1, conv([T(1), 1], [T(2), 1]));
    G.InputDelay = obiekt.tau;
    
    % Dyskretyzacja
    G_z = c2d(G, obiekt.Tp, 'zoh');
    G_z.Variable = 'z^-1';
    
    % Symulacja wyjścia Y
    Y = zeros(size(Y_step));
    Y(1:7) = Y_step(1:7); % załadowanie początkowych wartości (warunki początkowe)
    for k = 8:length(U)
        Y(k) = - G_z.Denominator{1}(2)*Y(k-1) - G_z.Denominator{1}(3)*Y(k-2) ...
            + G_z.Numerator{1}(2)*U(1, k-6) + G_z.Numerator{1}(3)*U(1, k-7) ...
            + G_z.Numerator{1}(2)*U(2, k-1) + G_z.Numerator{1}(3)*U(2, k-2);
    end
    
    E = sum((Y - Y_step).^2);
end