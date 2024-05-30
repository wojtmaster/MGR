%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kaskadowa regulacja poziomu wody w zbiornikach %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

%% Charakterystyka statyczna h_2(F)
% gdzie F = F_1 + F_D
% Autor: Wojciech Rogalski
% Data: 23.11.2023r.

%% Dane
P = 540;
C = 0.85;
alfa_1 = 26;
alfa_2 = 20;
Tp = 50;
steps = 100;

%% Punkt pracy
F_1 = linspace(0, 180, 1000);
F_10 = 100;
F_D = 30;
h_1 = ((F_1+F_D)/alfa_1).^2;
h_2 = ((F_1+F_D)/alfa_2).^2;
h_10 = ((F_10+F_D)/alfa_1)^2;
h_20 = ((F_10+F_D)/alfa_2)^2;

%% Prezentacja wyników
figure;
hold on;
plot(F_1, h_1, 'r');
plot(F_1, h_2, 'b');
plot(F_10, h_10, 'ro');
plot(F_10, h_20, 'bo');
hold off;
xlabel('F_{1}');
ylabel('h_{1}, h_{2}');
legend('h_{1}(F)', 'h_{2}(F)', 'Location', 'northwest');
title("Charakterystyka statyczna");