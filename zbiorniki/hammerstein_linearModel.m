%% LINEAR HAMMERSTEIN MODEL
clear all;
obiekt = Obiekt();
obiekt.linearization(90, 30);

% Zakres sterowania
F1_min = -45;
F1_max = 45;
FD_min = -15;
FD_max = 15;
F1_center = linspace(F1_min, F1_max, 5); % Środek zbiorów
FD_center = linspace(FD_min, FD_max, 3); % Środek zbiorów

%Generacja danych sterujących i RK4
% U = linspace(-45, 45, obiekt.kk);
% Y = ((90+U + obiekt.F_D0) / obiekt.alpha_2).^2;
% t = (0:length(U)-1) * obiekt.Tp; % Czas w sekundach (próbkowanie = 20s)

U = [linspace(-45, 45, 100);
    linspace(-15, 15, 100)];
[F1_grid, FD_grid] = meshgrid(U(1,:), U(2,:));  % kombinacje
Y = ((obiekt.F_10+F1_grid + obiekt.F_D0+FD_grid) / obiekt.alpha_2).^2;

surf(F1_grid+obiekt.F_10, FD_grid+obiekt.F_D0, Y);

%% Tworzenie początkowego systemu rozmytego TS
close all;
fis = sugfis('Name', 'Hammerstein', 'Type', 'sugeno');
fis = addInput(fis, [F1_min F1_max], 'Name', 'F1');
fis = addInput(fis, [FD_min FD_max], 'Name', 'FD');
rules_number = 15;
sigma_F1 = 12;
sigma_FD = 8;

% Definiowanie funkcji przynależności (gaussmf)
for i = 1:length(F1_center) 
    fis = addMF(fis, 'F1', 'gaussmf', [sigma_F1(i), F1_center(i)]);
end
for i = 1:length(FD_center) 
    fis = addMF(fis, 'FD', 'gaussmf', [sigma_FD, FD_center(i)]);
end

% figure;
% plotmf(fis, 'input', 1);
% figure;
% plotmf(fis, 'input', 2);

% Definiowanie wyjścia i początkowych następników (a_i * u + b_i)
fis = addOutput(fis, [-27 45], 'Name', 'Y_fuzzy');

% % Początkowe współczynniki (a_i, b_i)
% a_param = ones(1,rules_number)*0.6;
% b_param = ones(1,rules_number)*0.6;
% c_param = ones(1,rules_number)*1.9;

a_param = [	0.3677	0.4958	0.7417	0.4476	0.6090	0.7706	0.5915	0.5585	0.7681	0.5531	0.5410	0.6824	0.5186	0.7147	0.6810];
b_param = [	0.7801	0.8602	1.0268	0.5732	0.7117	0.7748	0.7151	0.5053	0.6096	0.6733	0.6520	0.5722	0.4357	0.5604	0.8445];
c_param = [	1.4143	2.2959	1.9559	1.8199	0.4147	1.6517	0.6147	0.3596	-0.6468	1.8442	0.9703	1.6752	3.1648	2.5298	2.8351];

% Parametry dla Y = [9 81]
% a_param = [0.5273    0.5902    0.5807    0.6200    0.6104    0.6144    0.6142    0.6565    0.5928    0.5546    0.5713    0.6247    0.6257    0.6126    0.6533];
% b_param = [0.5899    0.4921    0.6374    0.5616    0.6330    0.5702    0.5963    0.6288    0.7326    0.4502    0.5968    0.6048    0.6458    0.5991    0.5918];
% c_param = [42.6117   40.2665   37.8510   39.2992   37.6266   35.6608   36.3308   35.0988   35.2774   34.2218   37.7286   38.0913   38.2750   40.4708   44.0314];

% Dodanie reguł TS w postaci liniowej
for i = 1:rules_number
    fis = addMF(fis, 'Y_fuzzy', 'linear', [a_param(i), b_param(i), c_param(i)]);
end

% Reguły Takagi-Sugeno: [inputMF, outputMF, weight]
ruleList = [];
iter = 1;
for i = 1:length(F1_center)
    for j = 1:length(FD_center)
        ruleList(end+1, :) = [i, j, iter, 1, 1];
        iter = iter + 1;
    end
end

% Dodanie reguł do systemu
fis = addRule(fis, ruleList);

% % % Optymalizacja przy pomocy fminsearch
% initial_params = [a_param, b_param c_param];  % Początkowe wartości a_param i b_param
% options = optimset('Display', 'iter', 'MaxFunEvals', 2000, 'MaxIter', 1000); % Opcje optymalizacji
% optimal_params = fminsearch(@(params) linearCoeff(params, U, Y, F1_center, FD_center, sigma_F1, sigma_FD, rules_number), initial_params, options);

% % Po optymalizacji
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

%%
% Dane wejściowe (100x100 siatka → przekształć na 10000x2)
input_data = [F1_grid(:)+obiekt.F_10, FD_grid(:)+obiekt.F_D0];  % u1, u2 (każdy punkt z siatki)
output_data = Y(:);                    % odpowiadające pH (1D wektor)

clear fis;
% Liczba klastrów (zbiorów rozmytych) do znalezienia
% num_clusters = 25;  % Można zmieniać

fis = genfis2(input_data, output_data, 0.5);  % 0.5 = promień klastra (możesz go regulować)
% Wyświetl powierzchnię wyjściową
% gensurf(fis)

%% Check fuzzy static
for i = 1:length(fis.Rules)
    fis.Outputs.MembershipFunctions(i).Parameters(1) = a_optimal(i);
    fis.Outputs.MembershipFunctions(i).Parameters(2) = b_optimal(i);
    fis.Outputs.MembershipFunctions(i).Parameters(3) = c_optimal(i);
end

% gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));
% 
% sigma_F1 = 12;
% sigma_FD = 8;
% 
% for k = 1:length(U)
%     for l = 1:length(U)
%         deg_u1 = [gaussmf_val(U(1,l), sigma_F1, F1_center(1)), gaussmf_val(U(1,l), sigma_F1, F1_center(2)), ...
%             gaussmf_val(U(1,l), sigma_F1, F1_center(3)), gaussmf_val(U(1,l), sigma_F1, F1_center(4)), ...
%             gaussmf_val(U(1,l), sigma_F1, F1_center(5)),];
%         deg_u2 = [gaussmf_val(U(2,k), sigma_FD, FD_center(1)), gaussmf_val(U(2,k), sigma_FD, FD_center(2)), gaussmf_val(U(2,k), sigma_FD, FD_center(3))];
% 
%         degrees_all{k, l}(1) = deg_u1(1) * deg_u2(1);
%         degrees_all{k, l}(2) = deg_u1(1) * deg_u2(2);
%         degrees_all{k, l}(3) = deg_u1(1) * deg_u2(3);
%         degrees_all{k, l}(4) = deg_u1(2) * deg_u2(1);
%         degrees_all{k, l}(5) = deg_u1(2) * deg_u2(2);
%         degrees_all{k, l}(6) = deg_u1(2) * deg_u2(3);
%         degrees_all{k, l}(7) = deg_u1(3) * deg_u2(1);
%         degrees_all{k, l}(8) = deg_u1(3) * deg_u2(2);
%         degrees_all{k, l}(9) = deg_u1(3) * deg_u2(3);
%         degrees_all{k, l}(10) = deg_u1(4) * deg_u2(1);
%         degrees_all{k, l}(11) = deg_u1(4) * deg_u2(2);
%         degrees_all{k, l}(12) = deg_u1(4) * deg_u2(3);
%         degrees_all{k, l}(13) = deg_u1(5) * deg_u2(1);
%         degrees_all{k, l}(14) = deg_u1(5) * deg_u2(2);
%         degrees_all{k, l}(15) = deg_u1(5) * deg_u2(3);
%     end
% end
% 
% Y_out = zeros(size(Y));
% E_out = 0;
% % Oblicz odpowiedzi systemu rozmytego dla wszystkich U
% for i = 1:length(U)
%     for j = 1:length(U)
%         w = 0;
%         output = 0;
%         for k = 1:15
%             output = output + degrees_all{j,i}(k)*(a_optimal(k)*U(1,j) + b_optimal(k)*U(2,i) + c_optimal(k));
%             w = w + degrees_all{j, i}(k);
%         end
%         Y_out(i,j) = output / w;
%     end
% end

Y_out = zeros(size(Y));
for i = 1:length(U)
    for j = 1:length(U)
        Y_out(i,j) = evalfis(fis, [U(1,j), U(2,i)]) + obiekt.h_20;
    end
end

surf(F1_grid, FD_grid, Y_out);

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
T0 = [50, 50];
% Optymalna stała czasowa T = 193.4125      36.78401

% Szukanie optymalnej stałej czasowej
options = optimset('Display', 'off', 'MaxFunEvals', 2000, 'MaxIter', 1000); % Opcje optymalizacji
T_opt = fminsearch(fun, T0, options);

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

a = G_z.Denominator{1}(2:end);
b = G_z.Numerator{1}(2:end);

% plot(t, Y_step, 'r', 'LineWidth', 2);
% hold on;
% plot(t, Y, 'g');
% grid on;

%% Losowość
U = [repelem((rand(1, obiekt.kk/400) * 90 - 45), 400)
     repelem((rand(1, obiekt.kk/500) * 30 - 15), 500)];
U(1,1:8) = 0;
U(2,1:3) = 0;
[Y_real, Y_lin] = obiekt.rk4(U, obiekt.kk);
t = (0:length(U)-1) * obiekt.Tp; % Czas w sekundach (próbkowanie = 20s)

% Y_fuzzy = evalfis(fis, [U(1,:)', U(2,:)']);
% Y_0 = evalfis(fis, [10 -10]);

Y = zeros(1,obiekt.kk);
Y_fuzzy = zeros(1, obiekt.kk);

for k = 8:obiekt.kk
    Y_fuzzy(k) = evalfis(fis, [U(1,k), U(2,k)]);
    if(U(1,k-7) ~= 0)
        K = (Y_fuzzy(k-7)) / (U(1,k-7) + U(2,k-2));
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

%% DMC-SL
% Alokacja pamięci
y_zad = [repelem((rand(1, obiekt.kk/400) * 40 - 20), 400)];
y_zad(1:100) = 0;
u = zeros(2, obiekt.kk);
u(2,:) = [repelem((rand(1, obiekt.kk/500) * 30 - 15), 500)];

%% DMC(N, Nu, D, D_disturbance, lambda)
[s, s_disturbance] = obiekt.linearization(90, 30);

dmc = DMC(50, 4, 150, 150, 0.25);
dmc.dynamic_matrix(s);
dmc.past_matrix(s);
dmc.matrix_disturbance(s_disturbance);

delta_up = zeros(1, dmc.D-1);
delta_uz = zeros(1, dmc.D_disturbance);
delta_uk = 0;
delta_u = zeros(1, obiekt.kk);
F_10 = obiekt.F_10;
F_D0 = obiekt.F_D0;
h_20 = obiekt.h_20;
delay = obiekt.delay;
u(1,:) = zeros(1, obiekt.kk);

%% Sterowanie DMC
y_mod = zeros(size(y_zad));
y_mod(1:delay+2) = y_zad(1:delay+2);
y_fuzzy = zeros(size(y_zad));
for k = delay+3:obiekt.kk
    
    % if (mod(k,300) == 0)
    %     [s, s_disturbance] = obiekt.linearization(F_10 + u(1,k-1), F_D0 + u(2,k-1));
    %     dmc.dynamic_matrix(s);
    %     dmc.past_matrix(s);
    %     dmc.matrix_disturbance(s_disturbance);
    %     [a, b] = obiekt.sopdt();
    % end

    y_fuzzy(k-1) = evalfis(fis, [u(1,k-1), u(2,k-1)]);

    if(u(1,k-(delay+2)) + u(2,k-2) ~= 0)
        gain = y_fuzzy(k-(delay+2)) / (u(1,k-(delay+2)) + u(2,k-2));
    else
        gain = 1;
    end

    y_mod(k) = - a*[y_mod(k-1:-1:k-2)]' + gain*b*[u(1, k-(delay+1):-1:k-(delay+2))]' + gain*b*[u(2, k-1:-1:k-2)]';

    % Ograniczenia wartości sygnału wyjściowego, tj. wysokości h_2
    y_mod(k) = min(y_mod(k), dmc.y_max);
    y_mod(k) = max(y_mod(k), -dmc.y_max);

    % Przepisanie sterowań do wektora przeszłych sterowań
    delta_up = [delta_uk, delta_up(1:end-1)];
    delta_uz = [u(2, k) - u(2, k-1) , delta_uz(1:end-1)];

    % Oblicznie uchybu    
    e = y_zad(k) - y_mod(k);

    % Obliczenie przyrostu sterowania dla chwili (i+1) w chwili i-tej
    delta_uk = dmc.ke * e - dmc.ku * delta_up' - dmc.kz * delta_uz';
    
    % Ograniczenie wartości przyrostu sterowania
    delta_uk = min(delta_uk, dmc.delta_uk_max);
    delta_uk = max(delta_uk, -dmc.delta_uk_max);
    
    % Obliczenie sterowania dla chwili (i+1) w chwili i-tej
    u(1, k) = u(1, k-1) + delta_uk; 
    
    delta_F10 = obiekt.F_10 - F_10;

    % Ograniczenie sterowania
    u(1, k) = min(u(1, k), dmc.u_max);
    u(1, k) = max(u(1, k), -dmc.u_max);

    delta_u(k) = delta_uk;
end
E_u = dmc.lambda .* sum(delta_u.^2)/obiekt.kk;
E_y = sum((y_zad - y_mod).^2)/obiekt.kk;
E =  E_u + E_y;

fprintf("Błąd DMC \t E_u = %.3f \n", E_u);
fprintf("Błąd DMC \t E_y: %.3f \n", E_y);
fprintf("Błąd DMC \t E = %.3f \n\n", E);

% y_fig = figure;
figure(y_fig);
hold on;
stairs(0:obiekt.kk-1, y_mod, 'LineWidth', 0.8);
stairs(0:obiekt.kk-1, y_zad, 'r-', 'LineWidth', 0.8);
hold off;
% legend('y', 'y_{zad}', 'Location', 'best');
title(sprintf('Sygnał wyjściowy y(k)'));
xlabel('k');
ylabel('y(k)');
grid on;

% u_fig = figure;
figure(u_fig);
hold on;
stairs(0:obiekt.kk-1, u(1,:), 'LineWidth', 0.8);
hold off;
% legend('u', 'Location', 'best');
title(sprintf('Sygnał sterujący u(k)'));
xlabel('k');
ylabel('u(k)');
grid on;

%% Funkcja do optymalizacji współczynników
function E_out = linearCoeff(params, U, Y, F1_center, FD_center, sigma_F1, sigma_FD, rules_number)
    gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));

    % Parametry do optymalizacji
    a_param = params(1:rules_number);
    b_param = params(rules_number+1:2*rules_number);
    c_param = params(2*rules_number+1:end);

    for k = 1:length(U)
        for l = 1:length(U)
            % deg_u1 = [gaussmf_val(U(1,l), sigma_F1, F1_center(1)), gaussmf_val(U(1,l), sigma_F1, F1_center(2)), ...
            %     gaussmf_val(U(1,l), sigma_F1, F1_center(3)), gaussmf_val(U(1,l), sigma_F1, F1_center(4)), ...
            %     gaussmf_val(U(1,l), sigma_F1, F1_center(5)),];
            % deg_u2 = [gaussmf_val(U(2,k), sigma_FD, FD_center(1)), gaussmf_val(U(2,k), sigma_FD, FD_center(2)), gaussmf_val(U(2,k), sigma_FD, FD_center(3))];
            % 
            % degrees_all{k, l}(1) = deg_u1(1) * deg_u2(1);
            % degrees_all{k, l}(2) = deg_u1(1) * deg_u2(2);
            % degrees_all{k, l}(3) = deg_u1(1) * deg_u2(3);
            % degrees_all{k, l}(4) = deg_u1(2) * deg_u2(1);
            % degrees_all{k, l}(5) = deg_u1(2) * deg_u2(2);
            % degrees_all{k, l}(6) = deg_u1(2) * deg_u2(3);
            % degrees_all{k, l}(7) = deg_u1(3) * deg_u2(1);
            % degrees_all{k, l}(8) = deg_u1(3) * deg_u2(2);
            % degrees_all{k, l}(9) = deg_u1(3) * deg_u2(3);
            % degrees_all{k, l}(10) = deg_u1(4) * deg_u2(1);
            % degrees_all{k, l}(11) = deg_u1(4) * deg_u2(2);
            % degrees_all{k, l}(12) = deg_u1(4) * deg_u2(3);
            % degrees_all{k, l}(13) = deg_u1(5) * deg_u2(1);
            % degrees_all{k, l}(14) = deg_u1(5) * deg_u2(2);
            % degrees_all{k, l}(15) = deg_u1(5) * deg_u2(3);


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
    
    Y_out = zeros(size(Y));
    E_out = 0;
    % Oblicz odpowiedzi systemu rozmytego dla wszystkich U
    for i = 1:length(U)
        for j = 1:length(U)
            w = 0;
            output = 0;
            for k = 1:rules_number
                output = output + degrees_all{i,j}(k)*(a_param(k)*(U(1,j)) + b_param(k)*(U(2,i)) + c_param(k));
                w = w + degrees_all{i,j}(k);
            end
            Y_out(i,j) = output / w;
            % Y_out(i,j) = evalfis(fis, [U(1,j), U(2,i)]);
            % Funkcja kosztu - suma kwadratów błędów
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