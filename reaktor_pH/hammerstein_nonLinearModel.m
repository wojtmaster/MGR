%% LINEAR HAMMERSTEIN MODEL
clear;
obiekt = Obiekt();
%  pH_0, h_0, Q_10, Q_20, Q_30
obiekt.linearization(16.6, 0.55, 15.6);

%% Dane do wykresów 3D
% Zakresy sterowania
U_min = -15;
U_max = 15;
U = [linspace(U_min, U_max, 100);
    linspace(U_min, U_max, 100)];

[Q1_grid, Q3_grid] = meshgrid(U(1,:), U(2,:));  % kombinacje

% load pH_data.mat;
% load h_data.mat;

%% Przygotuj siatkę
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

figure;
surf(Q1_grid, Q3_grid, pH);
hold on;
xlabel('Q_1 [ml/s]');
ylabel('Q_3 [ml/s]');
zlabel('pH');
title('Wpływ dopływów Q_1 oraz Q_3 na stężenie substancji pH');
shading interp;
colorbar;

% q1_vals = linspace(min(Q1_grid(:)), max(Q1_grid(:)), 100);
% q3_vals = linspace(max(Q3_grid(:)), min(Q3_grid(:)), 100);
% 
% pH_cut = interp2(Q1_grid, Q3_grid, pH, q1_vals, q3_vals);
% 
% q1_new = linspace(min(Q1_grid(:)+5), max(Q1_grid(:)), 100);
% q3_new = linspace(max(Q3_grid(:)), min(Q3_grid(:)+5), 100);
% 
% pH_new = interp2(Q1_grid, Q3_grid, pH, q1_new, q3_new);
% 
% % Rysuj linię przekroju
% plot3(q1_vals, q3_vals, pH_cut, 'r');
% plot3(q1_new, q3_new, pH_new, 'g');

% saveas(gcf, 'D:/EiTI/MGR/raporty/raport_MGR/pictures/ph-_staticCharacteristic.png');  % Zapisuje jako plik PNG

figure;
surf(Q1_grid, Q3_grid, h);
xlabel('Q_1 [ml/s]');
ylabel('Q_3 [ml/s]');
zlabel('h [cm]');
title('Wpływ dopływów Q_1 oraz Q_3 na wysokość słupa cieczy h');
shading interp;
colorbar;
% view(45, 30);
% % saveas(gcf, 'D:/EiTI/MGR/raporty/raport_MGR/pictures/h_staticCharacteristic.png');  % Zapisuje jako plik PNG

%% Tworzenie początkowego systemu rozmytego TS
close all;
fis = sugfis('Name', 'Hammerstein', 'Type', 'sugeno');
rules_number = 8;
U_center = linspace(U_min, U_max, rules_number);
sigma = ones(1, rules_number) * 2.5;

fis = addInput(fis, [U_min U_max], 'Name', 'Q1');
fis = addInput(fis, [U_min U_max], 'Name', 'Q3');

% Definiowanie funkcji przynależności (gaussmf)
for i = 1:length(U_center) 
    fis = addMF(fis, 'Q1', 'gaussmf', [sigma(i), U_center(i)]);
    fis = addMF(fis, 'Q3', 'gaussmf', [sigma(i), U_center(i)]);
end

%%plotmf(fis, 'input', 1);

% Definiowanie wyjścia i początkowych następników (a_i * u + b_i)
fis = addOutput(fis, [U_min U_max], 'Name', 'Q_fuzzy');

% % Początkowe współczynniki (a_i, b_i, c_i)
% a_param = ones(1,rules_number)*-3.5;
% b_param = ones(1,rules_number)*0.2;
% % c_param = ones(1,rules_number)*1;

% % For h
% a_param = [	-1.4078	0.4264	1.4457	3.1947];
% b_param = [	-1.4918	1.3244	0.3397	3.1298];
% c_param = [	-26.0016	-5.4701	-6.3689	38.7067];

% For pH
a_param = [	-3.9 -3.5 -3.2 -3.5 0.1	-4.4 -3.7 -4.8];
b_param = [	0.01 0.4775	0.36 0.0454	0.1653	0.3265	0.3144	0.3512];

% % Optymalizacja przy pomocy fminsearch
% initial_params = [a_param, b_param];  % Początkowe wartości a_param i b_param
% options = optimset('Display', 'iter', 'MaxFunEvals', 2000, 'MaxIter', 1000); % Opcje optymalizacji
% optimal_params = fminsearch(@(params) nonlinearCoeff(params, U_center, U, pH, rules_number, sigma), initial_params, options);

% % Po optymalizacji
% a_optimal = optimal_params(1:rules_number);
% b_optimal = optimal_params(rules_number+1:2*rules_number);
% c_optimal = optimal_params(2*rules_number+1:end);
a_optimal = a_param;
b_optimal = b_param;
% % c_optimal = c_param;

% Wyświetlanie wyników optymalizacji
fprintf('a_param = [');
fprintf("\t%.4f", a_optimal);
fprintf('];\n');
fprintf('b_param = [');
fprintf("\t%.4f", b_optimal);
fprintf('];\n');
% fprintf('c_param = [');
% fprintf("\t%.4f", c_optimal);
% fprintf('];\n');

%% Check fuzzy static
clear deg;
gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));
% deg_u1 = zeros(1, sqrt(rules_number));
% deg_u2 = zeros(1, sqrt(rules_number));

clear degrees_all;
for k = 1:length(U)
    for l = 1:length(U)
        % for i = 1:sqrt(rules_number)
        %     deg_u1(i) = gaussmf_val(U(1,l), sigma(i), U_center(i));
        %     deg_u2(i) = gaussmf_val(U(2,k), sigma(i), U_center(i));
        % end
        % 
        % iter = 1;
        % for i = 1:sqrt(rules_number)
        %     for j = 1:sqrt(rules_number)
        %         degrees_all{k, l}(iter) = deg_u1(i) * deg_u2(j);
        %         iter = iter + 1;
        %     end
        % end

        for i = 1:rules_number
            % deg{k, l}(i) = gaussmf_val((U(1,l)-U(2,k)), sigma(i), U_center(i));
            deg{k, l}(i) = gaussmf_val((U(2,k)-U(1,l)), sigma(i), U_center(i));
        end
    end
end

Y_out = zeros(size(h));
for i = 1:length(U)
    for j = 1:length(U)
        output = 0;
        w = 0;
        for k = 1:rules_number
            % output = output + degrees_all{i,j}(k)*(a_optimal(k)*sinh(U(1,j)/7.5) + b_optimal(k)*sinh(U(2,i)/7.5) + c_optimal(k));
            % output = output + degrees_all{i,j}(k)*(a_optimal(k)*tanh(U(1,j)) + b_optimal(k)*tanh(U(2,i)) + c_optimal(k));
            % output = output + deg{i,j}(k)*(a_optimal(k)*tanh((U(1,j) - U(2,i))/2) + b_optimal(k));
            output = output + deg{i,j}(k)*(a_optimal(k)*tanh((U(2,i) - U(1,j))/2) + b_optimal(k));
            w = w + deg{i,j}(k);
        end
        Y_out(i,j) = output / w;
    end
end

E = sum(sum(pH - Y_out).^2 / obiekt.kk);
disp(E);

% Płaszczyzna 1 – Y_out – kolor czerwony
s1 = surf(Q1_grid, Q3_grid, Y_out);
hold on;
set(s1, 'FaceColor', 'red', 'EdgeColor', 'none'); % lub np. 'interp' dla interpolacji

% % Płaszczyzna 2 – pH – kolor niebieski
% s2 = surf(Q1_grid, Q3_grid, pH);
% set(s2, 'FaceColor', 'blue', 'EdgeColor', 'none');

xlabel('Q_1 [ml/s]');
ylabel('Q_3 [ml/s]');
zlabel('pH');
title('Wpływ dopływów Q_1 oraz Q_3 na stężenie substancji pH');
colorbar;

%% GA
rules_number = 8;
nvars = 2 * rules_number;
U_center = linspace(U_min, U_max, rules_number);
sigma = ones(1, rules_number) * 2.5;

% % % Zakresy graniczne (opcjonalnie, dopasuj do zakresu działania modelu)
lb = [ones(1, rules_number)*(-10)
      ones(1, rules_number)*(-10)];

ub = [ones(1, rules_number)*10
      ones(1, rules_number)*10];

% lb = [ones(1, rules_number)*(-3)
%       ones(1, rules_number)*(-3)
%       ones(1, rules_number)*(-10)];
% 
% ub = [ones(1, rules_number)*3
%       ones(1, rules_number)*3
%       ones(1, rules_number)*10];

% Definicja funkcji celu (loss function)
fitnessFcn = @(params) nonlinearCoeff(params, U_center, U, pH, rules_number, sigma);

% Opcje algorytmu genetycznego
opts = optimoptions('ga', ...
    'Display', 'iter', ...
    'MaxGenerations', 1000, ...
    'PopulationSize', nvars*10, ...
    'UseParallel', true, ...
    'PlotFcn', {@gaplotbestf}); % Włącza wykres błędu

% Wywołanie algorytmu genetycznego
[optimal_params, fval] = ga(fitnessFcn, nvars, [], [], [], [], lb, ub, [], opts);

% Po optymalizacji
a_optimal = optimal_params(1:rules_number);
b_optimal = optimal_params(rules_number+1:2*rules_number);
c_optimal = optimal_params(2*rules_number+1:end);

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

%% Check fuzzy static 1D
gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));

clear degrees_all;
for k = 1:length(U)
    l = length(U) - (k-1);

    % deg_u1 = zeros(1, sqrt(rules_number));
    % deg_u2 = zeros(1, sqrt(rules_number));
    for i = 1:rules_number
        deg{k}(i) = gaussmf_val((U(1,k)-U(2,l)), sigma(i), U_center(i));
        % deg_u1(i) = gaussmf_val(U(1,k), sigma(i), U_center(i));
        % deg_u2(i) = gaussmf_val(U(2,l), sigma(i), U_center(i));
    end

    % iter = 1;
    % for i = 1:sqrt(rules_number)
    %     for j = 1:sqrt(rules_number)
    %         degrees_all{k}(iter) = deg_u1(i) * deg_u2(j);
    %         iter = iter + 1;
    %     end
    % end
end

Y_out = zeros(size(pH_cut));
for i = 1:length(U)
    j = length(U) - (i-1);
    w = 0;
    output = 0;
    for k = 1:rules_number
        % index = round((k-1)/(rules_number-1) * (10-1) + 1);
        output = output + deg{i}(k)*(a_optimal(k)*tanh((U(1,i)- U(2,j))/2) + b_optimal(k));
        % output = output + degrees_all{i}(k)*(a_optimal(k)*tanh(U(1,i)) + b_optimal(k) * tanh(U(2,j)) + c_optimal(k));
        % output = output + degrees_all{i}(k)*(a_optimal(k)*U(1,i)) + b_optimal(k)*U(2,j) + c_optimal(k);
        w = w + deg{i}(k);
    end
    Y_out(i) = output / w;
end

figure;
plot(pH_new, 'r');
hold on;
plot(Y_out, 'b');



%% Identyfikacja dynamiki
[a_h, b_h, s_h] = obiekt.mse('h');
[a_pH, b_pH, s_pH] = obiekt.mse('pH');

%% Symulacja modelu Hammersteina
clear deg;
U = [repelem((rand(1, obiekt.kk/300) * 30 - 15), 300);
    zeros(1, obiekt.kk);
    repelem((rand(1, obiekt.kk/700) * 30 - 15), 700)];
[Y_real, Y_lin] = obiekt.modifiedEuler(U, obiekt.kk);
index = 2;
t = (0:length(U)-1) * obiekt.Tp; % Czas w sekundach (próbkowanie = 20s)

Y_out = zeros(1,obiekt.kk);
Y_fuzzy = zeros(1,obiekt.kk);
y_1 = zeros(1,obiekt.kk);
y_2 = zeros(1,obiekt.kk);
y_3 = zeros(1,obiekt.kk);

for k = 1:length(U)
    % deg_u1 = [gaussmf_val(U(1,k), sigma(1), U_center(1)), gaussmf_val(U(1,k), sigma(2), U_center(2)), ...
    %     gaussmf_val(U(1,k), sigma(3), U_center(3)), gaussmf_val(U(1,k), sigma(4), U_center(4)), ...
    %     gaussmf_val(U(1,k), sigma(5), U_center(5))];
    % deg_u2 = [gaussmf_val(U(3,k), sigma(1), U_center(1)), gaussmf_val(U(3,k), sigma(2), U_center(2)), ...
    %     gaussmf_val(U(3,k), sigma(3), U_center(3)), gaussmf_val(U(3,k), sigma(4), U_center(4)), ...
    %     gaussmf_val(U(3,k), sigma(5), U_center(5))];

    % iter = 1;
    % for i = 1:sqrt(rules_number)
    %     for j = 1:sqrt(rules_number)
    %         degrees_all(iter) = deg_u1(i) * deg_u2(j);
    %         iter = iter + 1;
    %     end
    % end

    for i = 1:rules_number
        deg(i) = gaussmf_val((U(3,k) - U(1,k)), sigma(i), U_center(i));            
    end

    output = 0;
    w = 0;
    for i = 1:rules_number
        output = output + deg(i)*(a_optimal(i)*tanh((U(3,k) - U(1,k))/2) + b_optimal(i));
        w = w + deg(i);
    end
    Y_fuzzy(k) = output / w;
end

for k = 4:obiekt.kk
    if(U(1,k-1) + U(3,k-1) ~= 0)
        K = Y_fuzzy(k-1) / (U(1,k-1) + U(3,k-1));
    else
        K = 1;
    end
    K_Q1 = Y_fuzzy(k-1) / (2*U(1,k-1));  % tylko jeśli U(1,k-1) ≠ 0
    K_Q3 = Y_fuzzy(k-1) / (2*U(3,k-1));

    % Y_out(k) = - a_h.Q1*Y_out(k-1) + K*b_h.Q1*U(1, k-1) + K*b_h.Q2*U(2, k-1) + K*b_h.Q3*U(3, k-1);

    y_1(k) = -a_pH.Q1 * y_1(k-1:-1:k-2)' + K*b_pH.Q1*U(1,k-1);
    y_2(k) = -a_pH.Q2 * y_2(k-1:-1:k-2)' + b_pH.Q2*U(2,k-1);
    y_3(k) = -a_pH.Q3 * y_3(k-1:-1:k-3)' + -K*b_pH.Q3*U(3,k-1);
    Y_out(k) = y_1(k) + y_2(k) + y_3(k);
end

% Wizualizacja wyników
figure;
plot(t, Y_real(index,:), 'b', t, Y_lin(index,:), 'g', t, Y_out, 'r');
legend('Euler', 'Euler liniowy', 'Hammerstein (optymalny TS)', 'Location', 'southwest');
title('Porównanie wyjścia układu rzeczywistego i modelu');
grid on;

figure;
subplot(1,3,1);
plot(t, y_1);
subplot(1,3,2);
plot(t, y_3);
subplot(1,3,3);
plot(t, Y_out);

E_lin = sum((Y_real(index,:) - Y_lin(index,:)).^2) / obiekt.kk;
E_out = sum((Y_real(index,:) - Y_out).^2) / obiekt.kk;
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
function E_out = nonlinearCoeff(params, U_center, U, Y, rules_number, sigma)
    gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));

    % Parametry do optymalizacji
    a_param = params(1:rules_number);
    b_param = params(rules_number+1:2*rules_number);
    % c_param = params(2*rules_number+1:end);

    % for k = 1:length(U)
    %     for l = 1:length(U)
    %         deg_u1 = zeros(1, sqrt(rules_number));
    %         deg_u2 = zeros(1, sqrt(rules_number));
    %         for i = 1:sqrt(rules_number)
    %             deg_u1(i) = gaussmf_val(U(1,l), sigma(i), U_center(i));
    %             deg_u2(i) = gaussmf_val(U(2,k), sigma(i), U_center(i));
    %         end
    % 
    %         iter = 1;
    %         for i = 1:sqrt(rules_number)
    %             for j = 1:sqrt(rules_number)
    %                 degrees_all{k, l}(iter) = deg_u1(i) * deg_u2(j);
    %                 iter = iter + 1;
    %             end
    %         end
    %     end
    % end

    % Y_out = zeros(size(Y));
    % E_out = 0;
    % % Oblicz odpowiedzi systemu rozmytego dla wszystkich U
    % for i = 1:length(U)
    %     for j = 1:length(U)
    %         w = 0;
    %         output = 0;
    %         for k = 1:rules_number
    %             % output = output + degrees_all{i,j}(k)*(a_param(k)*sinh(U(1,j)/7.5) + b_param(k)*sinh(U(2,i)/7.5) + c_param(k));
    %             output = output + degrees_all{i,j}(k)*(a_param(k)*tanh(U(1,j)) + b_param(k)*tanh(U(2,i)) + c_param(k));
    %             w = w + degrees_all{i,j}(k);
    %         end
    %         Y_out(i,j) = output / w;
    %         E_out = E_out + sum((Y(i,j) - Y_out(i,j))^2);
    %     end
    % end

    for k = 1:length(U)
        for l = 1:length(U)
            for i = 1:rules_number
                % deg{k,l}(i) = gaussmf_val((U(1,l)-U(2,k)), sigma(i), U_center(i));
                deg{k,l}(i) = gaussmf_val((U(2,k)-U(1,l)), sigma(i), U_center(i));
            end
        end
    end

    Y_out = zeros(size(Y));
    E_out = 0;
    % Oblicz odpowiedzi systemu rozmytego dla wszystkich U
    for i = 1:length(U)
        for j = 1:length(U)
            w = 0;
            output = 0;
            for k = 1:rules_number
                % output = output + degrees_all{i,j}(k)*(a_param(k)*sinh(U(1,j)/7.5) + b_param(k)*sinh(U(2,i)/7.5) + c_param(k));
                % output = output + deg{i,j}(k)*(a_param(k)*tanh((U(1,j)-U(2,i))/2) + b_param(k));
                output = output + deg{i,j}(k)*(a_param(k)*tanh((U(2,i)-U(1,j))/2) + b_param(k));
                w = w + deg{i,j}(k);
            end
            Y_out(i,j) = output / w;
            E_out = E_out + sum((Y(i,j) - Y_out(i,j))^2);
        end
    end
end