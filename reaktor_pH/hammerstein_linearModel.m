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

% Przygotuj siatkę
h = zeros(100,100);
pH = zeros(100,100);

% Petla po siatce sterowań
for i = 1:length(U)
    for j = 1:length(U)
        u = [ones(1, 200)*U(1,j);
             zeros(1, 200);
             ones(1, 200)*U(2,i)];
    
        [y, ~, ~] = obiekt.modifiedEuler(u, 200);
        h(i,j) =  y(1, end) + obiekt.h_0;
        pH(i,j) = y(2, end) + obiekt.pH_0; 
    end
end

% Rysuj 3D wykres
figure;
surf(Q1_grid, Q3_grid, pH);
xlabel('Q_3 [ml/s]');
ylabel('Q_1 [ml/s]');
zlabel('pH');
title('Wpływ dopływów Q_1 oraz Q_3 na stężenie substancji pH');
shading interp;
colorbar;
view(135, 0);
% saveas(gcf, 'D:/EiTI/MGR/raporty/raport_MGR/pictures/ph-_staticCharacteristic.png');  % Zapisuje jako plik PNG

figure;
surf(Q1_grid, Q3_grid, h);
xlabel('Q_1 [ml/s]');
ylabel('Q_3 [ml/s]');
zlabel('h [cm]');
title('Wpływ dopływów Q_1 oraz Q_3 na wysokość słupa cieczy h');
shading interp;
colorbar;
view(45, 30);
% saveas(gcf, 'D:/EiTI/MGR/raporty/raport_MGR/pictures/h_staticCharacteristic.png');  % Zapisuje jako plik PNG

%% 2. Tworzenie początkowego systemu rozmytego TS
close all;
fis = sugfis('Name', 'Hammerstein', 'Type', 'sugeno');
rules_number = 16;
U_center = linspace(U_min, U_max, sqrt(rules_number));
sigma = [4 4 4 4];

fis = addInput(fis, [U_min U_max], 'Name', 'Q1');
fis = addInput(fis, [U_min U_max], 'Name', 'Q3');

% Definiowanie funkcji przynależności (gaussmf)
for i = 1:length(U_center) 
    fis = addMF(fis, 'Q1', 'gaussmf', [sigma(i), U_center(i)]);
    fis = addMF(fis, 'Q3', 'gaussmf', [sigma(i), U_center(i)]);
end

% plotmf(fis, 'input', 1);

% Definiowanie wyjścia i początkowych następników (a_i * u + b_i)
fis = addOutput(fis, [U_min U_max], 'Name', 'Q_fuzzy');

% Początkowe współczynniki (a_i, b_i, c_i)
a_param = ones(1,rules_number)*0.001;
b_param = ones(1,rules_number)*0.25;
c_param = ones(1,rules_number)*0.2;

% % For h
% a_param = [0.3628    0.7838    1.0918    1.0118    0.9505    1.2950    1.1743    1.1477    1.4561];
% b_param = [0.7287    0.9401    1.0221    0.6728    1.1252    1.2062    1.1122    1.2291    1.6351];
% c_param = [21.5680   15.5153   15.0995   17.2411   14.0751   10.3066   11.4651   12.9631    6.7148];

% % For pH
% a_param = [0.0076	0.0111	0.0146	0.0071	0.0084	0.0051	0.0104	0.0085	0.0149	0.0184	0.0096	0.0031	0.0091	0.0131	0.0099	0.0057];
% b_param = [0.4562	0.1628	0.2621	0.2383	0.2939	0.3986	0.2130	0.1915	0.2072	0.3652	0.2975	0.2622	0.4104	0.3818	0.3222	0.2677];
% c_param = [10.9998	12.9682	9.1010	8.9379	7.7478	9.8031	12.0307	10.8059	6.1809	3.6940	4.4275	4.7971	6.3172	5.5293	0.7273	5.0499];

% Dodanie reguł TS w postaci liniowej
for i = 1:rules_number
    fis = addMF(fis, 'Q_fuzzy', 'linear', [a_param(i), b_param(i), c_param(i)]);
end

% Reguły Takagi-Sugeno: [inputMF, outputMF, weight]
ruleList = [];
iter = 1;
for i = 1:sqrt(rules_number)
    for j = 1:sqrt(rules_number)
        ruleList(end+1, :) = [i, j, iter, 1, 1];
        iter = iter + 1;
    end
end

fis = addRule(fis, ruleList);

% % Optymalizacja przy pomocy fminsearch
initial_params = [a_param, b_param c_param];  % Początkowe wartości a_param i b_param
options = optimset('Display', 'iter', 'MaxFunEvals', 4000, 'MaxIter', 2000); % Opcje optymalizacji
optimal_params = fminsearch(@(params) linearCoeff(params, U_center, U, pH, rules_number, sigma, obiekt), initial_params, options);

% Po optymalizacji
a_optimal = optimal_params(1:rules_number);
b_optimal = optimal_params(rules_number+1:2*rules_number);
c_optimal = optimal_params(2*rules_number+1:end);
% a_optimal = a_param;
% b_optimal = b_param;
% c_optimal = c_param;

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
input_data = [Q1_grid(:), Q3_grid(:)];  % u1, u2 (każdy punkt z siatki)
output_data = pH(:);                    % odpowiadające pH (1D wektor)

clear fis;
% Liczba klastrów (zbiorów rozmytych) do znalezienia
num_clusters = 150;  % Można zmieniać

fis = genfis2(input_data, output_data, 0.1);  % 0.5 = promień klastra (możesz go regulować)
% Wyświetl powierzchnię wyjściową
gensurf(fis)

%% Check fuzzy static
% for i = 1:rules_number
%     fis.Outputs.MembershipFunctions(i).Parameters(1) = a_optimal(i);
%     fis.Outputs.MembershipFunctions(i).Parameters(2) = b_optimal(i);
%     fis.Outputs.MembershipFunctions(i).Parameters(3) = c_optimal(i);
% end

Y_out = zeros(size(h));
for i = 1:length(U)
    for j = 1:length(U)
        Y_out(i,j) = evalfis(fis, [U(1,j), U(2,i)]);
    end
end

figure;
surf(Q1_grid, Q3_grid, Y_out);
% hold on;
% figure;
% surf(Q1_grid, Q3_grid, h);
% surf(Q1_grid, Q3_grid, pH);

%% Identyfikacja dynamiki
U = [ones(1, 100)
     zeros(1, 100)
     zeros(1, 100)];
[Y_step, ~] = obiekt.modifiedEuler(U, 100);
% Normalizacja odpowiedzi skokowej
Y_step(2,:) = Y_step(2,:) / Y_step(2,end);
t = (0:length(U)-1) * obiekt.Tp; % Czas w sekundach (próbkowanie = 20s)

% Funkcja celu (error)
fun = @(T) model_error(T, U(2,:), Y_step(2,:), obiekt);

% Startowa wartość (70)
T0 = 100;
% % For h
% T_Q1 = 180.5457;  % Q1
% T_Q2 = 180.5457;  % 
% T_Q3 = 180.5457;

% % For pH
T = 86.795;   % Q1
% T = 75.9905;  % Q2
% T = 143.4417;   % Q3 --> tu lepiej SOPDT

% % Szukanie optymalnej stałej czasowej
% options = optimset('Display', 'iter', 'MaxFunEvals', 2000, 'MaxIter', 1000); % Opcje optymalizacji
% T_opt = fminsearch(fun, T0, options);

% % Wyświetlenie wyniku
% disp(['Optymalna stała czasowa T = ', num2str(T_opt)]);

% G = tf(1, [T_opt, 1]);
G = tf(1, [T, 1]);
% G.InputDelay = obiekt.tau;
G_z = c2d(G, obiekt.Tp, 'zoh');
G_z.Variable = 'z^-1';

Y = zeros(1,100);
Y(1:7) = Y_step(2,1:7);
for k = 8:length(U)
    Y(k) = - G_z.Denominator{1}(2)*Y(k-1) + G_z.Numerator{1}(2)*U(1, k-1);
end

a.Q1 = G_z.Denominator{1}(2);
b.Q1 = G_z.Numerator{1}(2);
g_z.Q1 = G_z;

% plot(t, Y_step(2,:), 'r', 'LineWidth', 2);
% hold on;
% plot(t, Y, 'g');
% grid on;

%%
G1 = g_z.Q1;
G2 = g_z.Q2;
G3 = g_z.Q3;

% 1. Wyodrębnij liczniki i mianowniki
[b1, a1] = tfdata(G1, 'v');
[b2, a2] = tfdata(G2, 'v');
[b3, a3] = tfdata(G3, 'v');

% 2. Wyznacz wspólny mianownik jako najniższy wspólny wielokrotność (LCM) wielomianów
a_common = conv(a1, a3);  % może być później uproszczony

% 3. Przeskaluj liczniki tak, by odpowiadały wspólnemu mianownikowi
b1_common = conv(b1, a3);
b2_common = conv(b2, conv(a1, a3));
b3_common = conv(b3, a1);

% 4. Utwórz nowe transmitancje z tym samym mianownikiem
G1_new = tf(b1_common, a_common, G1.Ts);
G2_new = tf(b2_common, a_common, G2.Ts);
G3_new = tf(b3_common, a_common, G3.Ts);
G1_new.Variable = 'z^-1';
G2_new.Variable = 'z^-1';
G3_new.Variable = 'z^-1';

G1_new
G2_new
G3_new

A = G1_new.Denominator{1}(2:end);
B.Q1 = G1_new.Numerator{1}(2:end);
B.Q2 = G2_new.Numerator{1}(2:end);
B.Q3 = G3_new.Numerator{1}(2:end);

%% Symulacja modelu Hammersteina
% U = [repelem((rand(1, obiekt.kk/300) * 30 - 15), 300);
%     repelem((rand(1, obiekt.kk/700) * 10 - 5), 700);
%     repelem((rand(1, obiekt.kk/300) * 30 - 15), 300)];
U = [repelem((rand(1, obiekt.kk/300) * 30 - 15), 300);
    zeros(1, obiekt.kk);
    repelem((rand(1, obiekt.kk/300) * 30 - 15), 300)];
[Y_real, Y_lin] = obiekt.modifiedEuler(U, obiekt.kk);
t = (0:length(U)-1) * obiekt.Tp; % Czas w sekundach (próbkowanie = 20s)

Y_out = zeros(1,obiekt.kk);
Y_fuzzy = zeros(1,obiekt.kk);
y_1 = zeros(1,obiekt.kk);
y_2 = zeros(1,obiekt.kk);
y_3 = zeros(1,obiekt.kk);

Y_0 = evalfis(fis, [0 0]);

for k = 3:obiekt.kk
    Y_fuzzy(k) = evalfis(fis, [U(1,k), U(3,k)]);
    if(U(1,k-2) + U(3,k-2) ~= 0)
        K = (Y_fuzzy(k-2) - Y_0) / (U(1,k-2) + U(3,k-2));
    else
        K = 1;
    end
    % Y_out(k) = - a*Y_out(k-1) + K*b*U(1, k-1) + K*b*U(2, k-1) + K*b*U(3, k-1);
    % Y_out(k) = - A*Y_out(k-1:-1:k-3)' + K*B.Q1*U(1, k-1:-1:k-3)' + K*B.Q2*U(2, k-1:-1:k-3)' + K*B.Q3*U(3, k-1:-1:k-3)';
    % y_1(k) = -a.Q1 * y_1(k-1) + K*b.Q1*U(1,k-1);
    % y_2(k) = -a.Q2 * y_2(k-1) + K*b.Q2*U(2,k-1);
    % y_3(k) = -a.Q3 * y_3(k-1) + K*b.Q3*U(3,k-1);
    % Y_out(k) = y_1(k) + y_2(k) + y_3(k);
    Y_out(k) = - A*Y_out(k-1:-1:k-2)' + K*B.Q1*U(1, k-1:-1:k-2)' +  K*B.Q3*U(3, k-1:-1:k-2)';
end

% Wizualizacja wyników
figure;
plot(t, Y_real(1,:), 'b', t, Y_lin(1,:), 'g', t, Y_out, 'r');
legend('Euler', 'Euler liniowy', 'Hammerstein (optymalny TS)', 'Location', 'southwest');
title('Porównanie wyjścia układu rzeczywistego i modelu');
grid on;

E_lin = sum((Y_real(1,:) - Y_lin(1,:)).^2);
E_out = sum((Y_real(1,:) - Y_out).^2);
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
function E_out = linearCoeff(params, U_center, U, Y, rules_number, sigma, obiekt)
    gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));

    % Parametry do optymalizacji
    a_param = params(1:rules_number);
    b_param = params(rules_number+1:2*rules_number);
    c_param = params(2*rules_number+1:end);

    for k = 1:length(U)
        for l = 1:length(U)
            % deg_u1 = [gaussmf_val(U(1,l), sigma, U_center(1)), gaussmf_val(U(1,l), sigma, U_center(2)), ...
            %     gaussmf_val(U(1,l), sigma, U_center(3))];
            % deg_u2 = [gaussmf_val(U(2,k), sigma, U_center(1)), gaussmf_val(U(2,k), sigma, U_center(2)), ...
            %     gaussmf_val(U(2,k), sigma, U_center(3))];
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
            
            deg_u1 = zeros(1, sqrt(rules_number));
            deg_u2 = zeros(1, sqrt(rules_number));
            for i = 1:sqrt(rules_number)
                deg_u1(i) = gaussmf_val(U(1,l), sigma(i), U_center(i));
                deg_u2(i) = gaussmf_val(U(2,k), sigma(i), U_center(i));
            end

            iter = 1;
            for i = 1:sqrt(rules_number)
                for j = 1:sqrt(rules_number)
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
                output = output + degrees_all{i,j}(k)*(a_param(k)*U(1,j) + b_param(k)*U(2,i) + c_param(k));
                w = w + degrees_all{i,j}(k);
            end
            Y_out(i,j) = output / w;
            E_out = E_out + sum((Y(i,j) - Y_out(i,j))^2);
        end
    end
end

function E = model_error(T, U, Y_step, obiekt)
    % Tworzenie nowej transmitancji z aktualnym T
    G = tf(1, [T, 1]);
    % G.InputDelay = obiekt.tau;
    
    % Dyskretyzacja
    G_z = c2d(G, obiekt.Tp, 'zoh');
    G_z.Variable = 'z^-1';
    
    % Symulacja wyjścia Y
    Y = zeros(size(Y_step));
    Y(1:7) = Y_step(1:7); % załadowanie początkowych wartości (warunki początkowe)
    for k = 8:length(U)
        Y(k) = - G_z.Denominator{1}(2)*Y(k-1) + G_z.Numerator{1}(2)*U(k-1);
    end
    
    E = sum((Y - Y_step).^2);
end