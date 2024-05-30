function [] = modified_Euler(a, c, alpha_1, alpha_2, F_10, F_D0, h_10, h_20, V_10, V_20, t, kk, tau, T, u)

%% Skrypt przedstawiający zmodyfikowaną metodę Eulera

%% Alokacja pamięci
F_1in = zeros(1, kk);
F_D = zeros(1,kk);
y_1 = zeros(1, kk);
y_1L = zeros(1, kk);
y_2 = zeros(1, kk);
y_2L = zeros(1, kk);

%% Warunki początkowe
V_1 = V_10;
V_2 = V_20;
h_1 = h_10;
h_2 = h_20;
V_1e = V_10;
V_2e = V_20;
h_1e = h_10;
h_2e = h_20;
V_1L = V_10;
V_2L = V_20;
h_1L = h_10;
h_2L = h_20;
V_1eL = V_10;
V_2eL = V_20;
h_1eL = h_10;
h_2eL = h_20;

y_1(1) = h_10;
y_1L(1) = h_10;
y_2(1) = h_20;
y_2L(1) = h_20;

%% Funkcje
fun_1 = @(F_1, F_D, h_1) F_1 + F_D - alpha_1 * sqrt(h_1);
fun_1L = @(F_1, F_D, h_1) F_1 + F_D - alpha_1 * sqrt(h_10) - alpha_1 / (2*sqrt(h_10)) * (h_1-h_10);
fun_2 = @(h_1, h_2) alpha_1 * sqrt(h_1) - alpha_2 * sqrt(h_2); 
fun_2L = @(h_1, h_2) alpha_1 * sqrt(h_10) - alpha_2 * sqrt(h_20) + alpha_1 / (2*sqrt(h_10)) * (h_1-h_10) - alpha_2 / (2*sqrt(h_20)) * (h_2-h_20);

%% Wymuszenia
F_1in(1) = F_10;
F_1in(2 : 0.25*kk) = F_10 + u(1, 0.25*kk);
F_1in(0.25*kk+1 : 0.5*kk) = F_10 + u(1, 0.5*kk);
F_1in(0.5*kk+1 : 0.75*kk) = F_10 + u(1, 0.75*kk);
F_1in(0.75*kk+1 : end) = F_10 + u(1, end);

F_D(1) = F_D0;
F_D(2 : 0.25*kk) = F_D0+u(2, 0.25*kk);
F_D(0.25*kk+1 : 0.5*kk) = F_D0+u(2, 0.5*kk);
F_D(0.5*kk+1 : 0.75*kk) = F_D0+u(2, 0.75*kk);
F_D(0.75*kk+1 : end) = F_D0+u(2, end);

%% Modified Euler 
for i = 2:kk
    if i <= tau/T     
        V_1 = V_1 + T * fun_1(F_10, F_D(i), h_1);
        V_2 = V_2 + T * fun_2(h_1, h_2);
        h_1 = V_1 / a;
        h_2 = sqrt(V_2 / c);

        V_1e = V_1e + 1/2 * T * (fun_1(F_10, F_D(i), h_1e) + fun_1(F_10, F_D(i), h_1));
        V_2e = V_2e + 1/2 * T * (fun_2(h_1e, h_2e) + fun_2(h_1, h_2));

        V_1L = V_1L + T * fun_1L(F_10, F_D(i), h_1L);
        V_2L = V_2L + T * fun_2L(h_1L, h_2L);
        h_1L= V_1L / a;
        h_2L = sqrt(V_2L / c);

        V_1eL = V_1eL + 1/2 * T * (fun_1L(F_10, F_D(i), h_1eL) + fun_1L(F_10, F_D(i), h_1L));
        V_2eL = V_2eL + 1/2 * T * (fun_2L(h_1eL, h_2eL) + fun_2L(h_1L, h_2L));
    else                                      
        V_1 = V_1 + T * fun_1(F_1in(i - tau/T), F_D(i), h_1);
        V_2 = V_2 + T * fun_2(h_1, h_2);
        h_1 = V_1 / a;
        h_2 = sqrt(V_2 / c);

        V_1e = V_1e + 1/2 * T * (fun_1(F_1in(i - tau/T), F_D(i), h_1e) + fun_1(F_1in(i - tau/T), F_D(i), h_1));
        V_2e = V_2e + 1/2 * T * (fun_2(h_1e, h_2e) + fun_2(h_1, h_2));

        V_1L = V_1L + T * fun_1L(F_1in(i - tau/T), F_D(i), h_1L);
        V_2L = V_2L + T * fun_2L(h_1L, h_2L);
        h_1L = V_1L / a;
        h_2L = sqrt(V_2L / c);

        V_1eL = V_1eL + 1/2 * T * (fun_1L(F_1in(i - tau/T), F_D(i), h_1eL) + fun_1L(F_1in(i - tau/T), F_D(i), h_1L));
        V_2eL = V_2eL + 1/2 * T * (fun_2L(h_1eL, h_2eL) + fun_2L(h_1L, h_2L));
    end
        h_1e = V_1e / a;
        h_2e = sqrt(V_2e / c);
        h_1eL = V_1eL / a;
        h_2eL = sqrt(V_20 / c) + 1/(2*sqrt(V_20 * c)) * (V_2eL-V_20);

        y_1(i) = h_1e;
        y_1L(i) = h_1eL;
        y_2(i) = h_2e;
        y_2L(i) = h_2eL;
end

%% Prezentacja wyników
figure;
hold on;
plot(t, y_1, 'r','LineWidth',2);
plot(t, y_1L, 'b--','LineWidth',2);
hold off;
title('Wysokość słupa cieczy w zbiorniku 1. - h_1(t)');
legend('h_1', 'h_{1lin}');
xlabel('t [s]');
ylabel('h [cm]');
grid on;

figure;
hold on;
plot(t, y_2, 'g','LineWidth',2);
plot(t, y_2L, 'y--','LineWidth',2);
hold off;
title('Wysokość słupa cieczy w zbiorniku 2. - h_2(t)');
legend('h_2', 'h_{2lin}');
xlabel('t [s]');
ylabel('h [cm]');
grid on;
end