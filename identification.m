%% PPMGR
% Autor: Wojciech Rogalski
% Data: 11.03.2024r.
% Tytuł: Regulacja kaskadowa
clear all;
set(groot, 'DefaultTextFontName', 'Arial');
set(groot, 'DefaultAxesFontName', 'Arial');
set(groot, 'DefaultTextFontSize', 12);
set(groot, 'DefaultAxesFontSize', 10);

%% Dane
A = 540;
C = 0.85;
alpha_1 = 26;
alpha_2 = 20;

% Punkt pracy
F_10 = 90;
F_D0 = 30;
h_10 = ((F_10+F_D0)/alpha_1)^2;
h_20 = ((F_10+F_D0)/alpha_2)^2;
V_10 = A * h_10;
V_20 = C * h_20^2;

% Opóźnienie
tau = 100;
% Okres próbkowania
Tp = 20;
% Próbki dyskretne
kk = 5000;
% Wektor czasu
t = 0:Tp:(kk-1)*Tp;
% Zakres wyznaczanej charakterystyki statycznej
F_1_start = 45;
F_1_end = 135;
% Liczba podziałów wartości sterowania F_1 (do charakterystyki statycznej)
n = (F_1_end - F_1_start) * 2;

%% Charakrerystyka statyczna
[h_2, F_1] = static_characteristic(F_1_start, F_1_end, F_D0, alpha_2, n);

%% Podział danych statycznych
F_1_ucz = F_1(1:2:end-1);
F_1_wer = F_1(2:2:end);
h_2_ucz = h_2(1:2:end-1);
h_2_wer = h_2(2:2:end);

%% Fuzzy static charaacteristic - linear
[R, optimal_params] = static_characteristic_y_u(F_1, F_1_ucz, F_1_wer, h_2_ucz, h_2_wer, F_1_start, F_1_end, n, 'linear');

%% Fuzzy static charaacteristic - nonlinear
% [R, optimal_params] = static_characteristic_y_u(F_D0, alpha_2, F_1_start, F_1_end, n, 'nonlinear');

%% Wymuszenia + podział na zbiory danych dynamicznych
u = enforce(kk);
u_ucz = [u(1, 1:kk/2); u(2, 1:kk/2)];
u_wer = [u(1, kk/2+1:end); u(2, kk/2+1:end)];

%% Wyjście OBIEKTU - modified Euler + podział na zbiory danych dynamicznych
% Wyświetlane wartości na wykresach są sprowadzone do pkt. pracy --> wartości przyrostowe
y = modified_Euler(A, C, alpha_1, alpha_2, F_10, F_D0, h_20, V_10, V_20, t, kk, tau, Tp, u);
y_ucz = y(1:kk/2);
y_wer = y(kk/2+1:end);

%% Model o rzędzie dynamiki równym N
% Rząd dynamiki
N = 2;
% Opóźnienie
delay = tau/Tp;
Y_ucz = y_ucz(N+delay+1:end)';
M_ucz = u_ucz(1, N:end-(delay+1))';
M_wer = u_wer(1, N:end-(delay+1))';
% u(1, :)
for i = 2:N
    M_ucz = [M_ucz, u_ucz(1, (N+1)-i:end-(delay+i))'];
    M_wer = [M_wer, u_wer(1, (N+1)-i:end-(delay+i))'];
end
% u(2, :)
for i = 1:N
    M_ucz = [M_ucz, u_ucz(2, (N+delay+1)-i:end-i)'];
    M_wer = [M_wer, u_wer(2, (N+delay+1)-i:end-i)'];
end
% y(:)
for i = 1:N
    M_ucz = [M_ucz, y_ucz((N+delay+1)-i:end-i)'];
    M_wer = [M_wer, y_wer((N+delay+1)-i:end-i)'];
end

w = M_ucz\Y_ucz;

%% ARX
Y_ucz = [y_ucz(1:(N+delay))'; M_ucz*w];
Y_wer = [y_wer(1:(N+delay))'; M_wer*w];

E_ucz = sum((y_ucz'-Y_ucz).^2)/(kk/2);
E_wer = sum((y_wer' - Y_wer).^2)/(kk/2);
fprintf('Model ARX \n');
fprintf('E_ucz = %.3f \n', E_ucz);
fprintf('E_wer = %.3f \n', E_wer);

figure;
stairs(0:kk/2-1, Y_ucz, 'b', 'LineWidth', 1.2);
hold on;
stairs(0:kk/2-1, y_ucz, 'r', 'LineWidth', 0.8);
plot_title = sprintf('Zbiór uczący - y_{ucz}(k) \n N = %d     E_{ucz} = %.3f', N, E_ucz);
title(plot_title);
legend('y_{mod}', 'y');
xlabel('k');
ylabel('y(k)');
grid on;

figure;
stairs(0:kk/2-1, Y_wer, 'b', 'LineWidth', 1.2);
hold on;
stairs(0:kk/2-1, y_wer, 'r', 'LineWidth', 0.8);
plot_title = sprintf('Zbiór weryfikujący - y_{wer}(k) \n N = %d     E_{ucz} = %.3f', N, E_wer);
title(plot_title);
legend('y_{mod}', 'y');
xlabel('k');
ylabel('y(k)');
grid on;

%% OE
Y_ucz = zeros(kk/2, 1);
Y_wer = zeros(kk/2, 1);
Y_ucz(1:N+delay) = y_ucz(1:N+delay)';
Y_wer(1:N+delay) = y_wer(1:N+delay)';
for i = (N+delay)+1:kk/2
    M_temp_ucz = u_ucz(1, i-(delay+1));
    M_temp_wer = u_wer(1, i-(delay+1));
    for j = 2:N
        M_temp_ucz = [M_temp_ucz, u_ucz(1, i-(delay+j))];
        M_temp_wer = [M_temp_wer, u_wer(1, i-(delay+j))];
    end
    for j = 1:N
        M_temp_ucz = [M_temp_ucz, u_ucz(2, i-j)];
        M_temp_wer = [M_temp_wer, u_wer(2, i-j)];
    end
    for j = 1:N
        M_temp_ucz = [M_temp_ucz, Y_ucz(i-j)];
        M_temp_wer = [M_temp_wer, Y_wer(i-j)];
    end
    Y_ucz(i) = M_temp_ucz*w;
    Y_wer(i) = M_temp_wer*w;
end

E_ucz = sum((y_ucz'-Y_ucz).^2)/(kk/2);
E_wer = sum((y_wer' - Y_wer).^2)/(kk/2);
fprintf('Model OE \n');
fprintf('E_ucz = %.3f \n', E_ucz);
fprintf('E_wer = %.3f \n', E_wer);

figure;
stairs(0:kk/2-1, Y_ucz, 'b', 'LineWidth', 1.2);
hold on;
stairs(0:kk/2-1, y_ucz, 'r', 'LineWidth', 0.8);
plot_title = sprintf('Zbiór uczący - y_{ucz}(k) \n N = %d     E_{ucz} = %.3f', N, E_ucz);
title(plot_title);
legend('y_{mod}', 'y');
xlabel('k');
ylabel('y(k)');
grid on;

figure;
stairs(0:kk/2-1, Y_wer, 'b', 'LineWidth', 1.2);
hold on;
stairs(0:kk/2-1, y_wer, 'r', 'LineWidth', 0.8);
plot_title = sprintf('Zbiór weryfikujący - y_{wer}(k) \n N = %d     E_{ucz} = %.3f', N, E_wer);
title(plot_title);
legend('y_{mod}', 'y');
xlabel('k');
ylabel('y(k)');
grid on;

%%  Transmitancja
Tp = 20;
a = [1 -w(2*N+1:end)'];
b_1 = [0 w(1:N)'];
b_2 = [0 w(N+1:2*N)'];
G_z = tf({b_1, b_2}, a, 1, 'Variable', 'z', 'Ts', Tp);
G_z.Variable = 'z^-1';

a = a(2:end);
b(:,1) = b_1(2:end)/dcgain(G_z(1));
b(:,2) = b_2(2:end);
% b(:,2) = b_2(2:end)/dcgain(G_z(2));

y_mod_ucz = zeros(1, kk/2);
y_mod_wer = zeros(1, kk/2);
y_mod_ucz(1:N+delay) = y_ucz(1:N+delay);
y_mod_wer(1:N+delay) = y_wer(1:N+delay);

z = zeros(1, kk/2);

%% Model Hammersteina ARX
for k = N+delay+1:kk/2
    index = round((u_ucz(1,k)+F_10-F_1_start)*n / (F_1_end-F_1_start));
    z(k) = find_value(optimal_params, u_ucz(1,k)+F_10, index, R) - h_20;
    y_mod_ucz(k) = - a * [y_ucz(k-1:-1:k-N)]' + b(:, 1)'*[z(k-(delay+1):-1:k-(delay+N))]' + b(:,2)'*[u_ucz(2, k-1:-1:k-N)]';
end

% Zbiór weryfikujący
for k = N+delay+1:kk/2
    index = round((u_wer(1,k)+F_10-F_1_start)*n / (F_1_end-F_1_start));
    z(k) = find_value(optimal_params, u_wer(1,k)+F_10, index, R) - h_20;
    y_mod_wer(k) =  - a * [y_wer(k-1:-1:k-N)]' + b(:, 1)'*[z(k-(delay+1):-1:k-(delay+N))]' + b(:,2)'*[u_wer(2, k-1:-1:k-N)]';
end

E_ucz = sum((y_ucz - y_mod_ucz).^2)/(kk/2);
E_wer = sum((y_wer - y_mod_wer).^2)/(kk/2);
fprintf('Model ARX \n');
fprintf('E_ucz = %.3f \n', E_ucz);
fprintf('E_wer = %.3f \n', E_wer);

figure;
hold on;
stairs(0:kk/2-1, y_mod_ucz, 'b-', 'LineWidth', 1.2);
stairs(0:kk/2-1, y_ucz, 'r-', 'LineWidth', 0.8);
hold off;
xlabel('k');
ylabel('y_{ucz}(k)');
plot_title = sprintf('Zbiór uczący - y_{ucz}(k) \n E_{ucz} = %.3f', E_ucz);
title(plot_title);
legend('y_{mod}', 'y');
grid on;

figure;
hold on;
stairs(0:kk/2-1, y_mod_wer, 'b-', 'LineWidth', 1.2);
stairs(0:kk/2-1, y_wer, 'r-', 'LineWidth', 0.8);
hold off;
xlabel('k');
ylabel('y_{wer}(k)');
plot_title = sprintf('Zbiór weryfikujący - y_{wer}(k) \n E_{wer} = %.3f', E_wer);
title(plot_title);
legend('y_{mod}', 'y');
grid on;

%% Model Hammersteina OE
% Zbiór uczący
for k = N+delay+1:kk/2
    index = round((u_ucz(1,k)+F_10-F_1_start)*n / (F_1_end-F_1_start));
    z(k) = find_value(optimal_params, u_ucz(1,k)+F_10, index, R) - h_20;
    y_mod_ucz(k) = - a * [y_mod_ucz(k-1:-1:k-N)]' + b(:, 1)'*[z(k-(delay+1):-1:k-(delay+N))]' + b(:,2)'*[u_ucz(2, k-1:-1:k-N)]';
end

% Zbiór weryfikujący
for k = N+delay+1:kk/2
    index = round((u_wer(1,k)+F_10-F_1_start)*n / (F_1_end-F_1_start));
    z(k) = find_value(optimal_params, u_wer(1,k)+F_10, index, R) - h_20;
    y_mod_wer(k) =  - a * [y_mod_wer(k-1:-1:k-N)]' + b(:, 1)'*[z(k-(delay+1):-1:k-(delay+N))]' + b(:,2)'*[u_wer(2, k-1:-1:k-N)]';
end

E_ucz = sum((y_ucz - y_mod_ucz).^2)/(kk/2);
E_wer = sum((y_wer - y_mod_wer).^2)/(kk/2);
fprintf('Model OE \n');
fprintf('E_ucz = %.3f \n', E_ucz);
fprintf('E_wer = %.3f \n', E_wer);

figure;
hold on;
stairs(0:kk/2-1, y_mod_ucz, 'b-', 'LineWidth', 1.2);
stairs(0:kk/2-1, y_ucz, 'r-', 'LineWidth', 0.8);
hold off;
xlabel('k');
ylabel('y_{ucz}(k)');
plot_title = sprintf('Zbiór uczący - y_{ucz}(k) \n E_{ucz} = %.3f', E_ucz);
title(plot_title);
legend('y_{mod}', 'y');
grid on;

figure;
hold on;
stairs(0:kk/2-1, y_mod_wer, 'b-', 'LineWidth', 1.2);
stairs(0:kk/2-1, y_wer, 'r-', 'LineWidth', 0.8);
hold off;
xlabel('k');
ylabel('y_{wer}(k)');
plot_title = sprintf('Zbiór weryfikujący - y_{wer}(k) \n E_{wer} = %.3f', E_wer);
title(plot_title);
legend('y_{mod}', 'y');
grid on;

%% Model Wienera ARX
v = zeros(1, kk/2);

% Model ARX
% Zbiór uczący
for k = N+delay+1:kk/2
    v(k) = - a * [y_ucz(k-1:-1:k-N)]' + b(:, 1)'*[u_ucz(1, k-(delay+1):-1:k-(delay+N))]' + b(:,2)'*[u_ucz(2, k-1:-1:k-N)]';
    index = round((v(k)+F_10-F_1_start)*n / (F_1_end-F_1_start));
    y_mod_ucz(k) = find_value(optimal_params, v(k)+F_10, index, R) - h_20;
end

v = zeros(1, kk/2);
% Zbiór weryfikujący
for k = N+delay+1:kk/2
    v(k) =  - a * [y_wer(k-1:-1:k-N)]' + b(:, 1)'*[u_wer(1, k-(delay+1):-1:k-(delay+N))]' + b(:,2)'*[u_wer(2, k-1:-1:k-N)]';
    index = round((v(k)+F_10-F_1_start)*n / (F_1_end-F_1_start));
    y_mod_wer(k) = find_value(optimal_params, v(k)+F_10, index, R) - h_20;
end

E_ucz = sum((y_ucz - y_mod_ucz).^2)/(kk/2);
E_wer = sum((y_wer - y_mod_wer).^2)/(kk/2);
fprintf('Model ARX \n');
fprintf('E_ucz = %.3f \n', E_ucz);
fprintf('E_wer = %.3f \n', E_wer);

figure;
hold on;
stairs(0:kk/2-1, y_mod_ucz, 'b-', 'LineWidth', 1.2);
stairs(0:kk/2-1, y_ucz, 'r-', 'LineWidth', 0.8);
hold off;
xlabel('k');
ylabel('y_{ucz}(k)');
plot_title = sprintf('Zbiór uczący - y_{ucz}(k) \n E_{ucz} = %.3f', E_ucz);
title(plot_title);
legend('y_{mod}', 'y');
grid on;

figure;
hold on;
stairs(0:kk/2-1, y_mod_wer, 'b-', 'LineWidth', 1.2);
stairs(0:kk/2-1, y_wer, 'r-', 'LineWidth', 0.8);
hold off;
xlabel('k');
ylabel('y_{wer}(k)');
plot_title = sprintf('Zbiór weryfikujący - y_{wer}(k) \n E_{wer} = %.3f', E_wer);
title(plot_title);
legend('y_{mod}', 'y');
grid on;
%% Model WIenera OE
v = zeros(1, kk/2);

% Model ARX
% Zbiór uczący
for k = N+delay+1:kk/2
    v(k) = - a * [y_mod_ucz(k-1:-1:k-N)]' + b(:, 1)'*[u_ucz(1, k-(delay+1):-1:k-(delay+N))]' + b(:,2)'*[u_ucz(2, k-1:-1:k-N)]';
    index = round((v(k)+F_10-F_1_start)*n / (F_1_end-F_1_start));
    y_mod_ucz(k) = find_value(optimal_params, v(k)+F_10, index, R) - h_20;
end

v = zeros(1, kk/2);
% Zbiór weryfikujący
for k = N+delay+1:kk/2
    v(k) =  - a * [y_mod_wer(k-1:-1:k-N)]' + b(:, 1)'*[u_wer(1, k-(delay+1):-1:k-(delay+N))]' + b(:,2)'*[u_wer(2, k-1:-1:k-N)]';
    index = round((v(k)+F_10-F_1_start)*n / (F_1_end-F_1_start));
    y_mod_wer(k) = find_value(optimal_params, v(k)+F_10, index, R) - h_20;
end

E_ucz = sum((y_ucz - y_mod_ucz).^2)/(kk/2);
E_wer = sum((y_wer - y_mod_wer).^2)/(kk/2);
fprintf('Model ARX \n');
fprintf('E_ucz = %.3f \n', E_ucz);
fprintf('E_wer = %.3f \n', E_wer);

figure;
hold on;
stairs(0:kk/2-1, y_mod_ucz, 'b-', 'LineWidth', 1.2);
stairs(0:kk/2-1, y_ucz, 'r-', 'LineWidth', 0.8);
hold off;
xlabel('k');
ylabel('y_{ucz}(k)');
plot_title = sprintf('Zbiór uczący - y_{ucz}(k) \n E_{ucz} = %.3f', E_ucz);
title(plot_title);
legend('y_{mod}', 'y');
grid on;

figure;
hold on;
stairs(0:kk/2-1, y_mod_wer, 'b-', 'LineWidth', 1.2);
stairs(0:kk/2-1, y_wer, 'r-', 'LineWidth', 0.8);
hold off;
xlabel('k');
ylabel('y_{wer}(k)');
plot_title = sprintf('Zbiór weryfikujący - y_{wer}(k) \n E_{wer} = %.3f', E_wer);
title(plot_title);
legend('y_{mod}', 'y');
grid on;
