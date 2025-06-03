%% Valve types
u = linspace(0, 1, 1000);

% 1. Charakterystyka liniowa
y_linear = u;

% 2. Pierwiastek kwadratowy
y_sqrt = sqrt(u);

% 3. Equal percentage (równy procent) – model: y = (exp(k*u) - 1)/(exp(k) - 1)
k1 = 10;  % stała nieliniowości
y_eq_percent = (exp(k1 * u) - 1) / (exp(k1) - 1);

% 4. Charakterystyka hiperboliczna – model: y = u / (c + (1 - c)*u)
c = 5;
y_hyperbolic = u ./ (c + (1 - c) * u);

% Rysowanie wykresu
figure;
plot(u, y_linear, 'b', 'LineWidth', 2); hold on;
plot(u, y_sqrt, 'r', 'LineWidth', 2);
plot(u, y_eq_percent, 'g', 'LineWidth', 2);
plot(u, y_hyperbolic, 'm', 'LineWidth', 2);
grid on;

xlabel('Sygnał sterujący');
ylabel('Otwarcie zaworu');
title('Charakterystyki nieliniowe zaworów');
legend('Liniowa', 'Pierwiastek kw.', 'Równy procent', 'Hiperboliczna', 'Location', 'SouthEast');
axis([0 1 0 1]);

%% PT vs PTC vs NTC
% Zakres temperatur (w stopniach Celsjusza)
T_C = linspace(50, 150, 1000);  % Umiarkowany, realistyczny zakres
T_K = T_C + 273.15;            % Temperatura w Kelwinach (dla NTC)

% 1. PT100 (czujnik platynowy – liniowy model)
R0_PT = 100;                  % Rezystancja w 0°C
alpha_PT = 0.00385;           % Współczynnik temperaturowy
R_PT = R0_PT * (1 + alpha_PT * T_C);

% 2. PTC – nieliniowy wzrost (np. model wykładniczy)
R0_PTC = 100;
k_PTC = 0.02;                 % Stała wzrostu (dobrana do efektu)
R_PTC = R0_PTC * exp(k_PTC * T_C);

% 3. NTC – klasyczny model z równania Steinharta-Harta
R0_NTC = 10000;              % Rezystancja w 25°C
T0_K = 25 + 273.15;          % Referencyjna temperatura
B = 3950;                    % Stała B (typowa wartość dla NTC)
R_NTC = R0_NTC * exp(B * (1./T_K - 1/T0_K));

% Rysowanie wykresów
figure;
plot(T_C, R_PT, 'b', 'LineWidth', 2); hold on;
plot(T_C, R_PTC, 'r', 'LineWidth', 2);
plot(T_C, R_NTC, 'g', 'LineWidth', 2);
grid on;

xlabel('Temperatura [°C]');
ylabel('Rezystancja [\Omega]');
title('Charakterystyki statyczne PT, PTC i NTC');
legend('PT100', 'PTC', 'NTC', 'Location', 'NorthEast');

%% Hiperbolic fun
% Zakres x
x = linspace(-5, 5, 1000);

% Obliczanie funkcji hiperbolicznych
y_sinh = sinh(x);
y_cosh = cosh(x);
y_tanh = tanh(x);

% Rysowanie wykresów
figure;
plot(x, y_sinh, 'r', 'LineWidth', 2);
grid on;

xlabel('x');
ylabel('Wartość funkcji');
title('Funkcje hiperboliczna sinh(x)');
legend('sinh(x)', 'Location', 'NorthWest');

figure;
plot(x, y_cosh, 'b', 'LineWidth', 2);
grid on;

xlabel('x');
ylabel('Wartość funkcji');
title('Funkcje hiperboliczna cosh(x)');
legend('cosh(x)', 'Location', 'NorthWest');

figure;
plot(x, y_tanh, 'g', 'LineWidth', 2);
grid on;

xlabel('x');
ylabel('Wartość funkcji');
title('Funkcje hiperboliczna tanh(x)');
legend('tanh(x)', 'Location', 'NorthWest');
