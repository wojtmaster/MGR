clear;
% Przykładowe dane wejściowe i wyjściowe
u_start = 45;
u_end = 135;
u = linspace(u_start, u_end, (u_end-u_start)+1);
d = 30;
y = ((u+d) / 20).^2;

figure(1);
plot(u, y, 'b-');
hold on;

% Fuzzy sets
R = cell(1,5);
alpha = [-0.1, 0.1, 0.1, 0.1, 0.1];
c = [45 70 90 110 135];

for i = 1:length(R)
    for j = 1:length(u)
        if i == 1 || i == length(R)
            R{i}(j) = 1 / (1 + exp(-alpha(i)*(u(j)-c(i))));
        else
            R{i}(j) = 1 / (1 + exp(-alpha(i)*(u(j)-c(i-1)))) - 1 / (1 + exp(-alpha(i)*(u(j)-c(i+1))));
        end
    end
end

% figure;
% for i = 1:length(R)
%     plot(u, R{i});
%     hold on;
% end

a = ones(size(R));
b = ones(size(R));
 
% MNK: Minimalizacja sumy kwadratów błędów
fun = @(params) sum((y - fuzzy_linear_model(params, u, R)).^2);
initial_params = [a(1), b(1)];
for i = 2:length(R)
    initial_params = [initial_params, a(i), b(i)];
end
opt = optimset('MaxFunEvals', 10^6, 'MaxIter', 10^6);
optimal_params = fminsearch(fun, initial_params, opt);

% Prezentacja wyników
Y = fuzzy_linear_model(optimal_params, u, R);
figure(1);
plot(u, Y, 'ro');
axis([45 135 10 70]);
xlabel('u');
ylabel('y');
title('Charakterystyka statyczna y(u)');
legend('y(u)', 'y(u) - fuzzy', 'Location', 'northwest');

a = 540;
c = 0.85;
alpha_1 = 26;
alpha_2 = 20;

% Punkt pracy
F_10 = 90;
F_D0 = 30;
h_10 = ((F_10+F_D0)/alpha_1)^2;
h_20 = ((F_10+F_D0)/alpha_2)^2;
V_10 = a * h_10;
V_20 = c * h_20^2;

% Opóźnienie
tau = 100;
% Okres próbkowania
Tp = 20;
% Próbki dyskretne
kk = 500;
% Wektor czasu
t = 0:Tp:(kk-1)*Tp;

%% Wymuszenia
[u] = enforce(kk);

%% Transmitancje G(s) i G(z)
[G_s, G_z] = tf_function(a, c, alpha_1, alpha_2, V_10, V_20, tau, Tp);

%% Wartość zadana h_2
delta_h = set_value(kk, h_20);

%% Współczynniki do równań różnicowych na h_2
b(1) = [G_z.Numerator{2,1}(2)];
b(2) = [G_z.Numerator{2,1}(3)];

a(1) = [G_z.Denominator{2,1}(2)];
a(2) = [G_z.Denominator{2,1}(3)];

%% Regulator DMC
% Horyzonty
N = 100;
D = N;
Nu = 2;
lambda = 0.2;

% Ograniczenia (x_min = -x_max)
y_max = 10;
u_max = 15;
delta_u_max = 5;

% Odpwiedzi skokowe
s = step(G_z(2,1), Tp*(D-1));

%% Regulacja DMC - analitycznie
DMC_analitic(a, b, N, Nu, D, lambda, s, u, delta_h, y_max, u_max, delta_u_max, kk, tau/Tp, t);

% Znajdowanie wartości y na charakterystyce statycznej
u = 100;
y = find_value(optimal_params, u, u_start, R);
disp(y);