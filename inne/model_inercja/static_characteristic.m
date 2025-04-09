function [h_2, F_1] = static_characteristic(F_1_start, F_1_end, F_D0, alpha_2, n)

%% Charakterystyka statyczna h_2(F)
% gdzie F = F_1 + F_D
% Autor: Wojciech Rogalski
% Data: 23.11.2023r.

%% Punkt pracy
F_1 = linspace(F_1_start, F_1_end, n);
F_D = F_D0;
h_2 = ((F_1+F_D)/alpha_2).^2;

%% Prezentacja wyników
figure;
plot(F_1, h_2, 'b', 'LineWidth', 1.2);
xlabel('F_1');
ylabel('h_2');
xlim([45 135]);
legend('h_{2}(F_1)', 'Location', 'northwest');
title("Charakterystyka statyczna w zależności od F_1");
grid on;
end