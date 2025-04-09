%% MGR
% Autor: Wojciech Rogalski
% Data: 15.12.2024r.
% Tytuł: Porównanie modeli Hammersteina i Winera w regulacji rozmytej
clear all;
clc;
set(groot, 'DefaultTextFontName', 'Arial');
set(groot, 'DefaultAxesFontName', 'Arial');
set(groot, 'DefaultTextFontSize', 12);
set(groot, 'DefaultAxesFontSize', 10);

%% Inicjalizacja obiektu oraz linearyzacja
obiekt = Obiekt();

%% Testy 
% C_Ain_0, F_0, F_C_0, T_in_0, T_Cin_0, C_A0, C_B0, C_D0, T_0, T_C_0
obiekt.linearization(10, 0.1, 0.2, 300, 300, 4.8, 5.17, 0.03, 335.98, 321.86);
% C_Ain, F, T_in, F_C, T_Cin
u = [ones(1, obiekt.kk)*0
    ones(1, obiekt.kk)*0
    ones(1, obiekt.kk)*0
    ones(1, obiekt.kk)*-0.1
    ones(1, obiekt.kk)*0];
% y = obiekt.rk4(u);
y = obiekt.modifiedEuler(u);

if ~exist('y_figure', 'var') || ~isvalid(y_figure)
    y_figure = figure;
else
    figure(y_figure); % Jeśli istnieje, przełącz na nią
end

% Zrobić jeden wykres ze wszystkimi y
% Następnie zlinearyzować w tym "lepszym punkcie"
% Kilka wymuszeń dla każdego sterowania / wymuszenia

yyaxis left
h1 = plot(0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp, y(1,:), 'b-', 'LineWidth', 2);
ylabel('C');

yyaxis right
h2 = plot(0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp, y(2,:), 'r-', 'LineWidth', 2);
text(750, y(2, end)+0.05, ['F_C = ', sprintf('%.2f', u(4,1))]);
ylabel('T');

title('Wartości sygnałów wyjściowych y(k) - wymuszenie F_C');
xlabel('t [s]');
% Wymuszamy poprawną legendę, przypisując do odpowiednich uchwytów (h1, h2)
legend([h1, h2], {'C_D', 'T'}, 'Location', 'northwest');
hold on;
grid on;

% saveas(gcf, 'D:/EiTI/MGR/reaktor/wymuszenie_4.png');  % Zapisuje jako plik PNG

% obiekt.show_rk4(u, y);

%% Newton - Rapson
% Definicja danych (przykładowe wartości)
V = 10;        % m^3
V_C = 5;
rho = 1000;    % kg/m^3
rho_c = 1000;
c_p = 1;     % kJ/kg*K
c_pc = 4.2;
hA_C = 1300;   % W/m^2 K
T_C = 321.86;     % K
k_o1 = 4;      % 1/s
k_o2 = 172.2;  % 1/s
E_1 = 20900;   % J/mol
E_2 = 41800;   % J/mol
R = 8.3143;    % J/mol*K
H_1 = -41800;  % J/mol
H_2 = -83600;  % J/mol

% Funkcja do rozwiązywania układu równań
fun = @(x, u) [
    u(2)/V * (u(1) - x(1)) - k_o1 * x(1)^2 * exp(-E_1/(R*x(4)));
    -u(2)/V * x(2) + k_o1 * x(1)^2 * exp(-E_1/(R*x(4))) - k_o2 * x(2) * exp(-E_2/(R*x(4)));
    -u(2)/V * x(3) + k_o2 * x(2) * exp(-E_2/(R*x(4)));
    u(2)/V * (u(3) - x(4)) - 1/(rho * c_p) * (H_1 * k_o1 * x(1)^2 * exp(-E_1/(R*x(4))) + H_2 * k_o2 * x(2) * exp(-E_2/(R*x(4)))) - hA_C/(rho * c_p * V) * (x(4) - x(5));
    u(4)/V_C * (u(5) - x(5)) + hA_C/(rho_c *c_pc * V_C) * (x(4) - x(5));
];

% Początkowe przybliżenie
n = 100;
x0 = [4; 5; 0.01; 335; 321];
CA = zeros(1, n);
CB = zeros(1, n);
CD = zeros(1, n);
T = zeros(1, n);
TC = zeros(1, n);

CA_in = linspace(2, 18, n);
F = 0.1;
T_in = 300;
F_C = 0.2;
T_Cin = 300;

for i = 1:n
    % Rozwiązanie układu równań nieliniowych
    u = [CA_in(i), F(1), T_in(1), F_C(1), T_Cin(1)];

    fun_with_u = @(x) fun(x, u);

    % Rozwiązanie układu równań nieliniowych
    sol = fsolve(fun_with_u, x0);
    
    % Wyświetlenie wyników
    CA(i) = sol(1);
    CB(i) = sol(2);
    CD(i) = sol(3);
    T(i) = sol(4);
    TC(i) = sol(5);
end

subplot(2,1,1);
plot(CA_in, T, 'r-', 'LineWidth', 2);
title('T (C_{Ain})');
xlabel('C_{Ain}');
ylabel('T');
grid on;

subplot(2,1,2);
plot(CA_in, CD, 'b-', 'LineWidth', 2);
title('C_D (C_{Ain})');
xlabel('C_{Ain}');
ylabel('C_D');
grid on;

saveas(gcf, 'D:/EiTI/MGR/reaktor/static_char_1.png');  % Zapisuje jako plik PNG