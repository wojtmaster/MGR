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
% C_10, C_20, C_1in_0, C_2in_0, F_10, F_20, pH_10, pH_20, pH_0
obiekt.linearization(0.004, 0.004, 0.008, 0.008, 0.3, 0.3, 3.431, 11.903, 8.174);
% C_1in, F_1, C_2in, F_2
u = [ones(1, obiekt.kk)*0
    ones(1, obiekt.kk)*0
    ones(1, obiekt.kk)*0
    ones(1, obiekt.kk)*0];
% y = obiekt.rk4(u);
y = obiekt.modifiedEuler(u);

if ~exist('c_figure', 'var') || ~isvalid(c_figure)
    c_figure = figure;
else
    figure(c_figure); % Jeśli istnieje, przełącz na nią
end

plot(0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp, y(1,:), 'c-', 'LineWidth', 2);
hold on;
plot(0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp, y(2,:), 'm-', 'LineWidth', 2);
xlabel('t [s]');
ylabel('C');
title('Wartości sygnałów wyjściowych C_1, C_2 - wymuszenie F_2');
legend('C_1', 'C_2', 'Location', 'northwest');
grid on;
text(75, y(2, end), ['F_2 = ', sprintf('%.2f', u(4,1))]);

% saveas(gcf, 'D:/EiTI/MGR/reaktor_pH/wymuszenie_4_1.png');  % Zapisuje jako plik PNG

if ~exist('ph_figure', 'var') || ~isvalid(ph_figure)
    ph_figure = figure;
else
    figure(ph_figure); % Jeśli istnieje, przełącz na nią
end

plot(0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp, y(3,:), 'g-', 'LineWidth', 2);
text(75, y(3, end), ['F_2= ', sprintf('%.2f', u(4,1))]);
hold on;
xlabel('t [s]');
ylabel('pH');
title('Wartości sygnału wyjściowego pH - wymuszenie F_2');
legend('pH', 'Location', 'northwest');
grid on;

% saveas(gcf, 'D:/EiTI/MGR/reaktor_pH/wymuszenie_4_2.png');  % Zapisuje jako plik PNG

%% Charakterystyki statyczne
n = 100;
% C_1in, F_1, C_2in, F_2
u = [linspace(-0.008, 0.008, 100)
    linspace(-0.3, 0.3, 100)
    linspace(-0.008, 0.008, 100)
    linspace(-0.3, 0.3, 100)];
y = obiekt.staticCharakteristic(u, n);

%% Prezentacja wyników
figure;
plot(u(1,:), y(1,:), 'b-');
xlabel('C_{1in}');
ylabel('pH');
grid on;
title('Charakterystyka statyczna pH(C_{1in})');
% saveas(gcf, 'D:/EiTI/MGR/reaktor_pH/staticChar_1.png');  % Zapisuje jako plik PNG

figure;
plot(u(2,:), y(2,:), 'r-');
xlabel('F_1');
ylabel('pH');
grid on;
title('Charakterystyka statyczna pH(F_1)');
% saveas(gcf, 'D:/EiTI/MGR/reaktor_pH/staticChar_2.png');  % Zapisuje jako plik PNG

figure;
plot(u(3,:), y(3,:), 'b-');
hold on;
xlabel('C_{2in}');
ylabel('pH');
grid on;
title('Charakterystyka statyczna pH(C_{2in})');
% saveas(gcf, 'D:/EiTI/MGR/reaktor_pH/staticChar_3.png');  % Zapisuje jako plik PNG

figure;
plot(u(4,:), y(4,:), 'r-');
xlabel('F_2');
ylabel('pH');
grid on;
title('Charakterystyka statyczna pH(F_2)');
% saveas(gcf, 'D:/EiTI/MGR/reaktor_pH/staticChar_4.png');  % Zapisuje jako plik PNG