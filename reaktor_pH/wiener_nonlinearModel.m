%% LINEAR WIENER MODEL
clear;
obiekt = Obiekt();
%  pH_0, h_0, Q_10, Q_20, Q_30
obiekt.linearization(16.6, 0.55, 15.6);

[a_h, b_h, s_h] = obiekt.mse('h');
[a_pH, b_pH, s_pH] = obiekt.mse('pH');

%% Przygotowanie danych
U_min = -15;
U_max = 15;
U = [linspace(U_min, U_max, 100);
    linspace(U_min, U_max, 100)];
[Q1_grid, Q3_grid] = meshgrid(U(1,:), U(2,:));  % kombinacje

h = zeros(100,100);
pH = zeros(100,100);

% Petla po siatce sterowań
for i = 1:length(U)
    for j = 1:length(U)
        u = [ones(1, 200)*U(1,i);
             zeros(1, 200);
             ones(1, 200)*U(2,j)];
    
        [y, ~, ~] = obiekt.modifiedEuler(u, 200);
        h(i,j) =  y(1, end);
        pH(i,j) = y(2, end); 
    end
end

% Rysuj 3D wykres
figure;
surf(Q1_grid, Q3_grid, pH);
xlabel('Q_1 [ml/s]');
ylabel('Q_3 [ml/s]');
zlabel('pH');
title('Wpływ dopływów Q_1 oraz Q_3 na stężenie substancji pH');
% shading interp;
% colorbar;

figure;
surf(Q1_grid, Q3_grid, h);
xlabel('Q_1 [ml/s]');
ylabel('Q_3 [ml/s]');
zlabel('h [cm]');
title('Wpływ dopływów Q_1 oraz Q_3 na wysokość słupa cieczy h');
shading interp;
colorbar;

%% MSE
Y_h = zeros(100, 100);
Y_pH = zeros(100, 100);

y_h = zeros(1, 100);

y_1 = zeros(1, 100);
y_2 = zeros(1, 100);
y_pH = zeros(1, 100);

for i = 1:100
    for j = 1:100
        u = [ones(1,100) * U(1,i);
            ones(1,100) * U(2,j)];
        for k = 2:length(u)
            y_h(k) = - a_h.Q1*y_h(k-1) + b_h.Q1*u(1, k-1) + b_h.Q3*u(2, k-1);

            if k >= 4
                y_1(k) = -a_pH.Q1 * y_1(k-1:-1:k-2)' + b_pH.Q1*u(1,k-1);
                y_2(k) = -a_pH.Q3 * y_2(k-1:-1:k-3)' + b_pH.Q3*u(2,k-1);
                y_pH(k) = y_1(k) + y_2(k);
            end
        end
        Y_h(i,j) = y_h(end);
        Y_pH(i,j) = y_pH(end);
    end
end

figure;
surf(Q1_grid, Q3_grid, Y_h);
title('h');

figure;
surf(Q1_grid, Q3_grid, Y_pH);
title('pH');

%% Tworzenie systemu rozmytego TS dla Wienera
fis = sugfis('Name', 'Wiener', 'Type', 'sugeno');

Y_min = min(min(Y_h));
Y_max = max(max(Y_h));
rules_number = 7;
sigma = 5;
% Y_center = linspace(Y_min, Y_max, rules_number);
Y_center = [-19 -10 -4 0 4 10 19];
fis = addInput(fis, [Y_min Y_max], 'Name', 'y');

% Definiowanie funkcji przynależności (gaussmf)
for i = 1:length(Y_center) 
    fis = addMF(fis, 'y', 'gaussmf', [sigma, Y_center(i)]);
end

% plotmf(fis, 'input', 1);

% Definiowanie wyjścia i początkowych następników (a_i * y + b_i)
fis = addOutput(fis, [-5 5], 'Name', 'y_fuzzy');

% Początkowe współczynniki (a_i, b_i)
a_param = ones(1,rules_number)*10;
b_param = ones(1,rules_number)*2;

% % Początkowe współczynniki (a_i, b_i)
% a_param = ones(1,rules_number)*3.3;
% b_param = ones(1,rules_number)*0.1;

% % For h
% a_param = [	0.4405	18.6448	6.1357];
% b_param = [	-4.7708	8.3233	4.3980];

% For pH
a_param = [34.2411 -9.3827 0.1275 18.8864 -45.4857 152.1274	-31.5801];
b_param = [29.4027 -12.8669	1.9032 -15.3257	49.0425	-148.4382 35.4665];


% % Optymalizacja przy pomocy fminsearch
% initial_params = [a_param, b_param];  % Początkowe wartości a_param i b_param
% options = optimset('Display', 'iter', ...
%     'MaxFunEvals', 4000, ...
%     'MaxIter', 2000, ...
%     'TolFun', 1e-06); % Opcje optymalizacji
% optimal_params = fminsearch(@(params) nonlinearCoeff(params, Y_center, Y_h, h, sigma, rules_number), initial_params, options);

% % Po optymalizacji
% a_optimal = optimal_params(1:rules_number);
% b_optimal = optimal_params(rules_number+1:end);
a_optimal = a_param;
b_optimal = b_param;

% % Wyświetlanie wyników optymalizacji
% fprintf('a_param = [');
% fprintf("\t%.4f", a_optimal);
% fprintf('];\n');
% fprintf('b_param = [');
% fprintf("\t%.4f", b_optimal);
% fprintf('];\n');

% Check fuzzy static
gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));

Y_fuzzy = zeros(size(Y_h));
for i = 1:length(Y_h)
    for j = 1:length(Y_h)
        output = 0;
        w = 0;
        for k = 1:rules_number
            % degrees = gaussmf_val(Y_h(i,j), sigma, Y_center(k));
            % output = output + degrees*(a_optimal(k)*sinh(Y_h(i,j)/15) + b_optimal(k));
            degrees = gaussmf_val(Y_pH(i,j), sigma, Y_center(k));
            output = output + degrees*(a_optimal(k)*tanh(Y_pH(i,j)) + b_optimal(k));
            w = w + degrees;
        end

        Y_fuzzy(i,j) = output / w;
    end
end

E = sum(sum(pH - Y_fuzzy).^2) / obiekt.kk;
disp(E);

surf(Q1_grid, Q3_grid, Y_fuzzy);

%% GA
rules_number = 7;
nvars = 2 * rules_number;
% Y_center = linspace(Y_min, Y_max, rules_number);
Y_center = [-19 -10 -4 0 4 10 19];
sigma = 5;

% % % Zakresy graniczne (opcjonalnie, dopasuj do zakresu działania modelu)
lb = [ones(1, rules_number)*(-30)
      ones(1, rules_number)*(-30)];

ub = [ones(1, rules_number)*30
      ones(1, rules_number)*30];

% Definicja funkcji celu (loss function)
fitnessFcn = @(params) nonlinearCoeff(params, Y_center, Y_pH, pH, sigma, rules_number);

% Opcje algorytmu genetycznego
opts = optimoptions('ga', ...
    'Display', 'iter', ...
    'MaxGenerations', 2000, ...
    'PopulationSize', nvars*10, ...
    'UseParallel', true, ...
    'PlotFcn', {@gaplotbestf}); % Włącza wykres błędu

% Wywołanie algorytmu genetycznego
[optimal_params, fval] = ga(fitnessFcn, nvars, [], [], [], [], lb, ub, [], opts);

% Po optymalizacji
a_optimal = optimal_params(1:rules_number);
b_optimal = optimal_params(rules_number+1:2*rules_number);
% c_optimal = optimal_params(2*rules_number+1:end);

% Wyświetlanie wyników optymalizacji
fprintf('a_param = [');
fprintf("\t%.4f", a_optimal);
fprintf('];\n');
fprintf('b_param = [');
fprintf("\t%.4f", b_optimal);
fprintf('];\n');



%% Symulacja modelu Wienera
gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));

U = [repelem((rand(1, obiekt.kk/300) * 30 - 15), 300);
    zeros(1, obiekt.kk);
    repelem((rand(1, obiekt.kk/700) * 30 - 15), 700)];
U(1, 1:10) = 0;
U(3, 1:10) = 0;
t = (0:length(U)-1) * obiekt.Tp; % Czas w sekundach (próbkowanie = 20s)

[Y_real, Y_lin] = obiekt.modifiedEuler(U, obiekt.kk); % Symulacja rzeczywistego układu

Y_out = zeros(2, obiekt.kk);
Y_fuzzy = zeros(2, obiekt.kk);
y_1 = zeros(1,obiekt.kk);
y_2 = zeros(1,obiekt.kk);
y_3 = zeros(1,obiekt.kk);

index = 2;

for k = 4:obiekt.kk
    % Y_out(index, k) = - a_h.Q1*Y_out(index, k-1) + b_h.Q1*U(1,k-1) + b_h.Q2 * U(2, k-1) + b_h.Q3 * U(3,k-1);

    y_1(k) = -a_pH.Q1 * y_1(k-1:-1:k-2)' + b_pH.Q1*U(1,k-1);
    y_2(k) = -a_pH.Q2 * y_2(k-1:-1:k-2)' + b_pH.Q2*U(2,k-1);
    y_3(k) = -a_pH.Q3 * y_3(k-1:-1:k-3)' + b_pH.Q3*U(3,k-1);

    Y_out(index, k) = y_1(k) + y_2(k) + y_3(k);

    output = 0;
    w = 0;
    for i = 1:rules_number
        degrees = gaussmf_val(Y_out(index, k), sigma, Y_center(i));
        % output = output + degrees * (a_optimal(i)*sinh(Y_out(index, k)/15) + b_optimal(i));
        output = output + degrees * (a_optimal(i) * tanh(Y_out(index, k)) + b_optimal(i));
        w = w + degrees;
    end
    Y_fuzzy(index, k) = output / w;
end

% Wizualizacja wyników
figure;
plot(t, Y_real(index,:), 'b', t, Y_lin(index,:), 'g', t, Y_fuzzy(index, :), 'r');
legend('Euler', 'Euler liniowy', 'Wiener (optymalny TS)', 'Location', 'best');
title('Porównanie wyjścia układu rzeczywistego i modelu');
grid on;

E_lin = sum((Y_real(index,:) - Y_lin(index,:)).^2) / obiekt.kk;
E_out = sum((Y_real(index,:) - Y_fuzzy(index,:)).^2) / obiekt.kk;
fprintf("\nE_lin = %.3f\n", E_lin);
fprintf("E_out = %.3f\n", E_out);

%% Funkcja celu: oblicza błąd średniokwadratowy E_out
function E_out = nonlinearCoeff(params, U_center, U, Y, sigma, rules_number)
    gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));

    % Parametry do optymalizacji
    a_param = params(1:rules_number);  % a_param
    b_param = params(rules_number+1:end);  % b_param
    
    % Symulacja modelu Wienera
    Y_out = zeros(size(U));
    E_out = 0;

    for i = 1:length(U)
        for j = 1:length(U)
            w = 0;
            output = 0;
            for k = 1:rules_number
                degrees = gaussmf_val(U(i,j), sigma, U_center(k));
                output = output + degrees*(a_param(k)*sinh(U(i,j)/15) + b_param(k));
                % output = output + degrees*(a_param(k)*tanh(U(i,j)) + b_param(k));
                w = w + degrees;
            end
            Y_out(i,j) = output / w;
            E_out = E_out + sum((Y(i,j) - Y_out(i,j))^2);
        end
    end
end