function [] = static_characteristic(F_10, F_D0, h_10, h_20, alpha_1, alpha_2)

%% Charakterystyka statyczna h_2(F)
% gdzie F = F_1 + F_D
% Autor: Wojciech Rogalski
% Data: 23.11.2023r.

%% Punkt pracy
F_1 = linspace(0, 180, 1000);
F_D = F_D0;
h_1 = ((F_1+F_D)/alpha_1).^2;
h_2 = ((F_1+F_D)/alpha_2).^2;

%% Prezentacja wyników
figure;
subplot(2,1,1);
hold on;
plot(F_1, h_1, 'r');
plot(F_1, h_2, 'b');
plot(F_10, h_10, 'ro');
plot(F_10, h_20, 'bo');
hold off;
xlabel('F_{1}');
ylabel('h_{1}, h_{2}');
legend('h_{1}(F_1)', 'h_{2}(F_1)', 'Location', 'northwest');
title("Charakterystyka statyczna w zależności od F_1");
grid on;

%% Punkt pracy
F_1 = F_10;
F_D = linspace(0, 100, 1000);
h_1 = ((F_1+F_D)/alpha_1).^2;
h_2 = ((F_1+F_D)/alpha_2).^2;

%% Prezentacja wyników
subplot(2,1,2);
hold on;
plot(F_D, h_1, 'r');
plot(F_D, h_2, 'b');
plot(F_D0, h_10, 'ro');
plot(F_D0, h_20, 'bo');
hold off;
xlabel('F_{D}');
ylabel('h_{1}, h_{2}');
legend('h_{1}(F_D)', 'h_{2}(F_D)', 'Location', 'northwest');
title("Charakterystyka statyczna w zależności od F_D");
grid on;
end