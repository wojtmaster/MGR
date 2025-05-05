%% LINEAR HAMMERSTEIN MODEL
clear all;
obiekt = Obiekt();
obiekt.linearization(90, 30);

% Zakres sterowania
F1_min = -45;
F1_max = 45;
FD_min = -15;
FD_max = 15;
F1_center = linspace(F1_min, F1_max, 3); % Środek zbiorów
FD_center = linspace(FD_min, FD_max, 2); % Środek zbiorów

%Generacja danych sterujących i RK4
% U = linspace(-45, 45, obiekt.kk);
% Y = ((obiekt.F_10+U + obiekt.F_D0) / obiekt.alpha_2).^2;
% t = (0:length(U)-1) * obiekt.Tp; % Czas w sekundach (próbkowanie = 20s)

U = [linspace(-45, 45, 100);
    linspace(-15, 15, 100)];
[F1_grid, FD_grid] = meshgrid(U(1,:), U(2,:));  % kombinacje
Y = ((obiekt.F_10+F1_grid + obiekt.F_D0+FD_grid) / obiekt.alpha_2).^2 - obiekt.h_20;

%% Tworzenie początkowego systemu rozmytego TS
close all;
close all;
fis = sugfis('Name', 'Hammerstein', 'Type', 'sugeno');
fis = addInput(fis, [F1_min F1_max], 'Name', 'F1');
fis = addInput(fis, [FD_min FD_max], 'Name', 'FD');

% Definiowanie funkcji przynależności (gaussmf)
for i = 1:length(F1_center) 
    fis = addMF(fis, 'F1', 'gaussmf', [25, F1_center(i)]);
end
for i = 1:length(FD_center) 
    fis = addMF(fis, 'FD', 'gaussmf', [15, FD_center(i)]);
end
figure;
plotmf(fis, 'input', 1);
figure;
plotmf(fis, 'input', 2);

%% Tworzenie początkowego systemu rozmytego TS
sigma_F1 = 25;
sigma_FD = 15;
rules_number = 6;
% Początkowe współczynniki (a_i, b_i)
% a_param = ones(1,6)*9;
% b_param = ones(1,6)*3;
% c_param = ones(1,6)*1.9;

a_param = [	5.5084	5.2897	16.4147	18.9656	4.5280	7.4589];
b_param = [	0.6532	0.8794	0.6018	0.9663	2.3825	2.6390];
c_param = [	1.1386	8.8990	-11.4339	7.0581	3.7385	5.8944];

% a_param = [3.1360    2.2899   10.7762   10.6849    4.6489    5.9102];
% b_param = [0.3367    0.6262    1.2488    1.8362    1.1683    1.3966];
% c_param = [24.0212   30.2372   29.2767   41.4704   39.3933   57.1661];

% % Optymalizacja przy pomocy fminsearch
% initial_params = [a_param, b_param c_param];  % Początkowe wartości a_param i b_param
% options = optimset('Display', 'iter', 'MaxFunEvals', 5000, 'MaxIter', 2500); % Opcje optymalizacji
% optimal_params = fminsearch(@(params) nonlinearCoeff(params, sigma_F1, sigma_FD, F1_center, FD_center, U, Y, rules_number), initial_params, options);

% Po optymalizacji
% a_optimal = optimal_params(1:rules_number);
% b_optimal = optimal_params(rules_number+1:2*rules_number);
% c_optimal = optimal_params(2*rules_number+1:end);
a_optimal = a_param;
b_optimal = b_param;
c_optimal = c_param;

% Wyświetlanie wyników optymalizacji
fprintf('a_param = [');
fprintf("\t%.4f", a_optimal);
fprintf('];\n');
fprintf('b_param = [');
fprintf("\t%.4f", b_optimal);
fprintf('];\n');
fprintf('c_param = [');
fprintf("\t%.4f", c_optimal);
fprintf('];\n');

%% Check fuzzy static
gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));

for k = 1:length(U)
    for l = 1:length(U)
        deg_u1 = [gaussmf_val(U(1,l), sigma_F1, F1_center(1)), gaussmf_val(U(1,l), sigma_F1, F1_center(2)), gaussmf_val(U(1,l), sigma_F1, F1_center(3))];
        deg_u2 = [gaussmf_val(U(2,k), sigma_FD, FD_center(1)), gaussmf_val(U(2,k), sigma_FD, FD_center(2))];

        degrees_all{k, l}(1) = deg_u1(1) * deg_u2(1);
        degrees_all{k, l}(2) = deg_u1(1) * deg_u2(2);
        degrees_all{k, l}(3) = deg_u1(2) * deg_u2(1);
        degrees_all{k, l}(4) = deg_u1(2) * deg_u2(2);
        degrees_all{k, l}(5) = deg_u1(3) * deg_u2(1);
        degrees_all{k, l}(6) = deg_u1(3) * deg_u2(2);
    end
end

Y_out = zeros(size(Y));
for i = 1:length(U)
    for j = 1:length(U)
        w = 0;
        output = 0;
        for k = 1:rules_number
            output = output + degrees_all{i,j}(k)*(a_optimal(k)*sinh(U(1,j)/22.5) + b_optimal(k)*sinh(U(2,i)/7.5) + c_optimal(k));
            w = w + degrees_all{i,j}(k);
        end
        Y_out(i,j) = output / w + obiekt.h_20;
    end
end

figure;
surf(F1_grid, FD_grid, Y_out);
% hold on;
% surf(F1_grid, FD_grid, Y);

%% Identyfikacja dynamiki
U = [ones(1, 100)
    zeros(1, 100)];
[Y_step, ~] = obiekt.rk4(U, 100);
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

% plot(t, Y_step, 'b', t, Y, 'g');
% grid on;

%% Losowość
U = [repelem((rand(1, obiekt.kk/400) * 90 - 45), 400)
     repelem((rand(1, obiekt.kk/500) * 30 - 15), 500)];
[Y_real, Y_lin] = obiekt.rk4(U, obiekt.kk);
t = (0:length(U)-1) * obiekt.Tp; % Czas w sekundach (próbkowanie = 20s)

Y = zeros(1,obiekt.kk);
Y_fuzzy = zeros(1,obiekt.kk);

clear degrees_all;
deg_u1 = [gaussmf_val(0, sigma_F1, F1_center(1)), gaussmf_val(0, sigma_F1, F1_center(2)), gaussmf_val(0, sigma_F1, F1_center(3))];
deg_u2 = [gaussmf_val(0, sigma_FD, FD_center(1)), gaussmf_val(0, sigma_FD, FD_center(2))];

degrees_all(1) = deg_u1(1) * deg_u2(1);
degrees_all(2) = deg_u1(1) * deg_u2(2);
degrees_all(3) = deg_u1(2) * deg_u2(1);
degrees_all(4) = deg_u1(2) * deg_u2(2);
degrees_all(5) = deg_u1(3) * deg_u2(1);
degrees_all(6) = deg_u1(3) * deg_u2(2);

w = 0;
output = 0;
for k = 1:rules_number
    output = output + degrees_all(k)*(a_optimal(k)*0 + b_optimal(k)*0 + c_optimal(k));
    w = w + degrees_all(k);
end
Y_0 = output / w;

for k = 1:obiekt.kk
    deg_u1 = [gaussmf_val(U(1,k), sigma_F1, F1_center(1)), gaussmf_val(U(1,k), sigma_F1, F1_center(2)), gaussmf_val(U(1,k), sigma_F1, F1_center(3))];
    deg_u2 = [gaussmf_val(U(2,k), sigma_FD, FD_center(1)), gaussmf_val(U(2,k), sigma_FD, FD_center(2))];
    
    degrees_all(1) = deg_u1(1) * deg_u2(1);
    degrees_all(2) = deg_u1(1) * deg_u2(2);
    degrees_all(3) = deg_u1(2) * deg_u2(1);
    degrees_all(4) = deg_u1(2) * deg_u2(2);
    degrees_all(5) = deg_u1(3) * deg_u2(1);
    degrees_all(6) = deg_u1(3) * deg_u2(2);
    
    w = 0;
    output = 0;
    for i = 1:rules_number
        output = output + degrees_all(i)*(a_optimal(i)*sinh(U(1,k)/22.5) + b_optimal(i)*sinh(U(2,k)/7.5) + c_optimal(i));
        w = w + degrees_all(i);
    end
    Y_fuzzy(k) = output / w;
end

for k = 8:obiekt.kk
    if(U(1,k-7) ~= 0)
        K = (Y_fuzzy(k-7)) / (U(1,k-7) + U(2, k-2));
    else
        K = 1;
    end
    Y(k) = - G_z.Denominator{1}(2)*Y(k-1) - G_z.Denominator{1}(3)*Y(k-2) ...
        + K * G_z.Numerator{1}(2)*U(1, k-6) + K * G_z.Numerator{1}(3)*U(1, k-7) ...
        + K * G_z.Numerator{1}(2)*U(2, k-1) + K * G_z.Numerator{1}(3)*U(2, k-2);
end

plot(t, Y_real, 'b');
hold on;
plot(t, Y_lin, 'r');
plot(t, Y, 'g');
grid on;
legend('Y_{real}', 'Y_{lin}', 'Y');

%% Funkcja do optymalizacji współczynników
function E_out = nonlinearCoeff(params, sigma_F1, sigma_FD, F1_center, FD_center, U, Y, rules_number)
    gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));

    % Parametry do optymalizacji
    a_param = params(1:rules_number);
    b_param = params(rules_number+1:2*rules_number);
    c_param = params(2*rules_number+1:end);

    for k = 1:length(U)
        for l = 1:length(U)
            deg_u1 = zeros(1, length(F1_center));
            deg_u2 = zeros(1, length(FD_center));
            for i = 1:length(F1_center)
                deg_u1(i) = gaussmf_val(U(1,l), sigma_F1, F1_center(i));
            end
            for i = 1:length(FD_center)
                deg_u2(i) = gaussmf_val(U(2,k), sigma_FD, FD_center(i));
            end

            iter = 1;
            for i = 1:length(F1_center)
                for j = 1:length(FD_center)
                    degrees_all{k, l}(iter) = deg_u1(i) * deg_u2(j);
                    iter = iter + 1;
                end
            end
        end
    end
    
    Y_out = zeros(size(U));
    E_out = 0;
    % Oblicz odpowiedzi systemu rozmytego dla wszystkich U
    for i = 1:length(U)
        for j = 1:length(U)
            w = 0;
            output = 0;
            for k = 1:rules_number
                output = output + degrees_all{i,j}(k)*(a_param(k)*sinh(U(1,j)/22.5) + b_param(k)*sinh(U(2,i)/7.5) + c_param(k));
                w = w + degrees_all{i,j}(k);
            end
            Y_out(i,j) = output / w;
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