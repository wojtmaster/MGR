%% Projekt
clear all;

%% Dane
P = 540;
C = 0.85;
alfa_1 = 26;
alfa_2 = 20;
Tp = 50;
steps = 100;

%% Punkt pracy
F_10 = 100;
F_D0 = 30;
h_10 = ((F_10+F_D0)/alfa_1)^2;
h_20 = ((F_10+F_D0)/alfa_2)^2;
V_10 = P * h_10;
V_20 = C * h_20^2;
tau = 100;

%% Alokacja pamięci
F_1 = zeros(1 ,steps);
F_1in = zeros(1, steps);
F_D = zeros(1, steps);
F_D(1:steps) = F_D0;
V_1F = zeros(1, steps);
V_2F = zeros(1, steps);
h_1F = zeros(1, steps);
h_2F = zeros(1, steps);
V_1F(1:tau/Tp+1) = V_10;
V_2F(1:tau/Tp+1) = V_20;
h_1F(1) = h_10;
h_2F(1) = h_20;
t = zeros(1, steps);

%% Rozmywanie F_1
sets = 5;
n = 1000;
width = 3;
F_1Fmax = 120;
F_1Fmin = 80;
F_1F = zeros(1, sets);
F_1F(1) = 80;
F_1F(2) = 90;
F_1F(3) = 100;
F_1F(4) = 110;
F_1F(5) = 120;
x = linspace(80, F_1Fmax, n);
w_arr = cell(1,5);

for i = 1:sets
    w_arr{i} = gaussmf(x, [width, F_1F(i)]);
end

figure;
hold on;
plot(x, w_arr{1});
plot(x, w_arr{2});
plot(x, w_arr{3});
plot(x, w_arr{4});
plot(x, w_arr{5});
hold off;
title('Funkcje przynależności wartości F_{1}');
xlabel('F_{1}');
ylabel('mi(F_{1})');
axis([80 120 0 1]);

%% Punkty pracy dla różnych wartości F_1F
V_1F0 = zeros(1, sets);
V_2F0 = zeros(1, sets);
h_1F0 = zeros(1, sets);
h_2F0 = zeros(1, sets);
V_1eLF = zeros(1, sets);
V_2eLF = zeros(1, sets);
h_1eLF = zeros(1, sets);
h_2eLF = zeros(1, sets);
w = zeros(1, sets);

for i = 1:sets
    h_1F0(i) = ((F_1F(i) + F_D0) / alfa_1)^2;
    h_2F0(i) = ((F_1F(i) + F_D0) / alfa_2)^2;
    V_1F0(i) = P * h_1F0(i);
    V_2F0(i) = C * (h_2F0(i))^2;
end

%% Funkcje różniczkujące
fun_1L = @(h_1, h10, F_1, F_D) F_1 + F_D - alfa_1 * sqrt(h10) - alfa_1 / (2*sqrt(h10)) * (h_1-h10);
fun_2L = @(h_1, h10, h_2, h20) alfa_1 * sqrt(h10) - alfa_2 * sqrt(h20) + alfa_1 / (2*sqrt(h10)) * (h_1-h10) - alfa_2 / (2*sqrt(h20)) * (h_2-h20);

F_1(1) = F_10;
F_1in(1:40) = F_10;
F_1in(41:steps) = F_10;

for i = 2:steps
    if(i < tau/Tp+1)
        F_1(i) = F_10;
        h_1F(i) = h_10;
        h_2F(i) = h_20;
        t(i) = (i-1)*Tp;
        continue;
    else
        F_1(i) = F_1in(i-tau/Tp);
    end
    number = round((F_1(i) - F_1Fmin) * (n - 0)/(F_1Fmax - F_1Fmin) + 0);
    for j = 1:sets
        w(j) = w_arr{j}(number);
        if w(j) < 1e-05
            V_1eLF(j) = 0;
            h_1eLF(j) = 0;
            V_2eLF(j) = 0;
            h_2eLF(j) = 0;
            continue;
        else
            V_1LF = V_1F(i-1) + Tp * fun_1L(h_1F(i-1), h_1F0(j), F_1(i-1), F_D(i-1));
            V_2LF = V_2F(i-1) + Tp * fun_2L(h_1F(i-1), h_1F0(j), h_2F(i-1), h_2F0(j));
            h_1LF = V_1LF / P;
            h_2LF = sqrt(V_2LF / C);
        
            V_1eLF(j) = V_1F(i-1) + 1/2 * Tp * (fun_1L(h_1F(i-1), h_1F0(j), F_1(i-1), F_D(i-1)) + fun_1L(h_1LF, h_1F0(j), F_1F(j), F_D(i)));
            V_2eLF(j) = V_2F(i-1) + 1/2 * Tp * (fun_2L(h_1F(i-1), h_1F0(j), h_2F(i-1), h_2F0(j)) + fun_2L(h_1LF, h_1F0(j), h_2LF, h_2F0(j)));
            h_1eLF(j) = V_1eLF(j) / P;
            h_2eLF(j) = sqrt(V_2eLF(j) / C);
        end
    end

        V_1_sum = 0;
        V_2_sum = 0;
        h_1_sum = 0;
        h_2_sum = 0;
        w_sum = 0;

    for j = 1:sets
        V_1_sum = V_1_sum + V_1eLF(j) * w(j);
        V_2_sum = V_2_sum + V_2eLF(j) * w(j);
        h_1_sum = h_1_sum + h_1eLF(j) * w(j);
        h_2_sum = h_2_sum + h_2eLF(j) * w(j);
        w_sum = w_sum + w(j);
    end

    V_1F(i) = V_1_sum / w_sum;
    V_2F(i) = V_2_sum / w_sum;
    h_1F(i) = h_1_sum / w_sum;
    h_2F(i) = h_2_sum / w_sum;
    t(i) = (i-1)*Tp;

end

%% Prezentacja wyników
figure;
subplot(2,1,1);
hold on;
plot(t, h_1F, 'r','LineWidth',2);
hold off;
title('h_{1}(t)');
legend('h_1');
xlabel('t [s]');
ylabel('h [cm]');
grid on;

subplot(2,1,2);
hold on;
plot(t, h_2F, 'g','LineWidth',2);
hold off;
title('h_{2}(t)');
legend('h_2');
xlabel('t [s]');
ylabel('h [cm]');
grid on;