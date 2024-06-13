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
kk = 150;
% Wektor czasu
t = 0:Tp:(kk-1)*Tp;

%% Wymuszenia + podział na zbiory danych dynamicznych
u = enforce(kk);
% close;
u_ucz = [u(1, 1:kk/2); u(2, 1:kk/2)];
u_wer = [u(1, kk/2+1:end); u(2, kk/2+1:end)];

%% Wyjście OBIEKTU - modified Euler + podział na zbiory danych dynamicznych
% Wyświetlane wartości na wykresach są sprowadzone do pkt. pracy --> wartości przyrostowe
y = modified_Euler(A, C, alpha_1, alpha_2, F_10, F_D0, h_20, V_10, V_20, t, kk, tau, Tp, u);

%% 
fig = figure;
figure(fig);
plot(t, y, 'r-','LineWidth',2);
hold on;
grid on;

%% Model inercjalny
T_0 = tau;
K_0 = 0.6025;
T_1 = 212;
T_2 = 15;
s = tf('s');
num = K_0;
den =conv([T_1 1], [T_2 1]);

G_s = tf(num, den, 'InputDelay', T_0);

G_z = c2d(G_s, Tp, 'zoh');
G_z.Variable = 'z^-1';


figure(fig);
step(G_z);

a(1:2) = G_z.Denominator{1}(2:end);
b(1:2) = G_z.Numerator{1}(2:end);