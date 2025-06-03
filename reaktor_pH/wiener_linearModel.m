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
             ones(1, 200)*0.5;
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
shading interp;
colorbar;
view(-45, 30);

figure;
surf(Q1_grid, Q3_grid, h);
xlabel('Q_1 [ml/s]');
ylabel('Q_3 [ml/s]');
zlabel('h [cm]');
title('Wpływ dopływów Q_1 oraz Q_3 na wysokość słupa cieczy h');
shading interp;
colorbar;

%% MSE - h
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

Y_min = min(min(Y_pH));
Y_max = max(max(Y_pH));
rules_number = 15;
sigma = 2;
% Y_center = linspace(-25, 25, rules_number);
Y_center = [-25 -15 -10 -8 -6 -4 -2 0 2 4 6 8 10 15 25];
fis = addInput(fis, [Y_min Y_max], 'Name', 'y');

% Definiowanie funkcji przynależności (gaussmf)
for i = 1:length(Y_center) 
    fis = addMF(fis, 'y', 'gaussmf', [sigma, Y_center(i)]);
end

% plotmf(fis, 'input', 1);

% Definiowanie wyjścia i początkowych następników (a_i * y + b_i)
fis = addOutput(fis, [-5 5], 'Name', 'y_fuzzy');

% Początkowe współczynniki (a_i, b_i)
a_param = ones(1,rules_number)*0.1;
b_param = ones(1,rules_number)*0.1;

% % For h
% a_param = [-0.0366 0.2409 0.5852 1.1168];
% b_param = [-15.5957 -4.4613	4.9938 7.2292];

% % For pH

% Dodanie reguł TS w postaci liniowej
for i = 1:length(Y_center)
    fis = addMF(fis, 'y_fuzzy', 'linear', [a_param(i), b_param(i)]);
end

% Reguły Takagi-Sugeno: [inputMF, outputMF, weight]
ruleList = [];
iter = 1;
for i = 1:rules_number
    ruleList(end+1, :) = [i, i, 1, 1];
    iter = iter + 1;
end

% Dodanie reguł do systemu
fis = addRule(fis, ruleList);

% Optymalizacja przy pomocy fminsearch
initial_params = [a_param, b_param];  % Początkowe wartości a_param i b_param
options = optimset('Display', 'final', 'MaxFunEvals', 4000, 'MaxIter', 2000); % Opcje optymalizacji
optimal_params = fminsearch(@(params) linearCoeff(params, Y_center, Y_pH, pH, sigma, rules_number), initial_params, options);

% Po optymalizacji
a_optimal = optimal_params(1:rules_number);
b_optimal = optimal_params(rules_number+1:end);
% a_optimal = a_param;
% b_optimal = b_param;

% Wyświetlanie wyników optymalizacji
fprintf('a_param = [');
fprintf("\t%.4f", a_optimal);
fprintf('];\n');
fprintf('b_param = [');
fprintf("\t%.4f", b_optimal);
fprintf('];\n');

%% ANFIS
% Dane wejściowe
% X = Y_pH(:);
X = tanh(Y_pH(:)/30);
Y = pH(:);
data_train = [X Y];

% Ustawienia
numMFs = 4;                 % liczba funkcji przynależności na wejście
mfType = 'gaussmf';         % typ MF: gaussmf, gbellmf, trimf, ...

% Generowanie systemu rozmytego
fis = genfis1(data_train, numMFs, mfType, 'linear');

% Trening ANFIS
options = anfisOptions('InitialFIS', fis, 'EpochNumber', 200, ...
    'DisplayANFISInformation', 0, 'DisplayErrorValues', 0);
fis_trained = anfis(data_train, options);

% Predykcja
Y_pred = evalfis(fis_trained, X);
pH_pred = reshape(Y_pred, size(Q1_grid));

% Rysunek
figure;
surf(Q1_grid, Q3_grid, pH_pred);
xlabel('Q_1'); ylabel('Q_3'); zlabel('pH');
title('Dokładniejsza charakterystyka statyczna TS (genfis1 + ANFIS)');
shading interp; colormap turbo;

%% GENFIS
input_data = Y_pH(:);
output_data = pH(:);

% opt = genfisOptions('GridPartition');   % lub 'SubtractiveClustering'
% opt.InputMembershipFunctionType = 'gbellmf';  % dobre do lokalnych zmian
% opt.NumMembershipFunctions = 11;         % lub więcej, jeśli masz dużo danych
% opt.OutputMembershipFunctionType = 'linear';  % TS system = Sugeno 1. rzędu

opt = genfisOptions('SubtractiveClustering');
opt.ClusterInfluenceRange = 0.1;  % im mniejsze, tym więcej reguł
opt.SquashFactor = 1.25;  % tłumienie wpływu sąsiednich klastrów
opt.AcceptRatio = 0.6;
opt.RejectRatio = 0.4;

fis = genfis(input_data, output_data, opt);

for i = 1:length(Y_pH)
    for j = 1:length(Y_pH)
        Y_out(i,j) = evalfis(fis, Y_pH(i,j)) + obiekt.pH_0;
    end
end

surf(Q1_grid, Q3_grid, Y_out);
xlabel('Q_1'), ylabel('Q_3'), zlabel('pH');

%% Check fuzzy static
for i = 1:length(fis.Rules)
    fis.Outputs.MembershipFunctions(i).Parameters(1) = a_optimal(i);
    fis.Outputs.MembershipFunctions(i).Parameters(2) = b_optimal(i);
end

Y_fuzzy = zeros(size(Y_h));
for i = 1:length(Y_h)
    for j = 1:length(Y_h)
        Y_fuzzy(i,j) = evalfis(fis, Y_h(i,j)) + obiekt.pH_0;
    end
end

surf(Q1_grid, Q3_grid, Y_fuzzy);

%% Symulacja modelu Wienera
gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));

% U = [repelem((rand(1, obiekt.kk/300) * 30 - 15), 300);
%     zeros(1, obiekt.kk);
%     repelem((rand(1, obiekt.kk/700) * 30 - 15), 700)];
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
    % Y_out(index, k) = - a_h.Q1*Y_out(k-1)' + b_h.Q1*U(1,k-1) + b_h.Q2 * U(2, k-1) + b_h.Q3 * U(3,k-1);
    % Y_fuzzy(index, k) = evalfis(fis, Y_out(index, k));
    
    y_1(k) = -a_pH.Q1 * y_1(k-1:-1:k-2)' + b_pH.Q1*U(1,k-1);
    y_2(k) = -a_pH.Q2 * y_2(k-1:-1:k-2)' + b_pH.Q2*U(2,k-1);
    y_3(k) = -a_pH.Q3 * y_3(k-1:-1:k-3)' + b_pH.Q3*U(3,k-1);

    Y_out(index, k) = y_1(k) + y_2(k) + y_3(k);
    Y_fuzzy(index, k) = evalfis(fis, Y_out(index, k));
end

% Wizualizacja wyników
figure;
plot(t, Y_real(index,:), 'b', t, Y_lin(index,:), 'g', t, Y_fuzzy(index, :), 'r');
legend('Euler', 'Euler liniowy', 'Wiener (optymalny TS)', 'Location', 'best');
title('Porównanie wyjścia układu rzeczywistego i modelu');
grid on;

E_lin = sum((Y_real(index,:) - Y_lin(index,:)).^2);
E_out = sum((Y_real(index,:) - Y_fuzzy(index,:)).^2);
fprintf("\nE_lin = %.3f\n", E_lin);
fprintf("E_out = %.3f\n", E_out);

%% Funkcja do optymalizacji współczynników
function E_out = linearCoeff(params, U_center, U, Y, sigma, rules_number)
    gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));

    % Parametry do optymalizacji
    a_param = params(1:rules_number);
    b_param = params(rules_number+1:end);

    Y_out = zeros(size(U));
    E_out = 0;

    for i = 1:length(U)
        for j = 1:length(U)
            for k = 1:rules_number
                degrees(k) = [gaussmf_val(U(i,j), sigma, U_center(k))];
            end

            w = 0;
            output = 0;
            for k = 1:rules_number
                output = output + degrees(k)*(a_param(k)*U(i,j) + b_param(k));
                w = w + degrees(k);
            end
            Y_out(i,j) = output / w;
            E_out = E_out + sum((Y(i,j) - Y_out(i,j))^2);
        end
    end
end