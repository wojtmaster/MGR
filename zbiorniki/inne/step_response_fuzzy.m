%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Regulacja przepływu wody między dwoma zbiornikami %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Skrypt generujący odpowiedzi skokowe dla wymuszenia jednostkowego
clear all;

%% Dane
P = 540;
C = 0.85;
alfa_1 = 26;
alfa_2 = 20;
Tp = 50;    
k_end = 200;

%% Punkt pracy
F_10 = 90;
F_D0 = 30;
h_10 = ((F_10+F_D0)/alfa_1)^2;
h_20 = ((F_10+F_D0)/alfa_2)^2;
V_10 = P * h_10;
V_20 = C * h_20^2;
tau = 100;

%% Równania stanu
A = [-alfa_1/(2*P*sqrt(h_10)) 0; alfa_1/(4*C*h_20*sqrt(h_10)) -alfa_2/(4*C*h_20*sqrt(h_20))];
B = [1/P 1/P; 0 0];
C = [0 1];

s1 = tf('s');
G_s = C*(s1*eye(2) - A)^(-1)*B;
G_z = c2d(G_s, Tp, 'tustin');

%% DMC parameters
N = 210;    %horyzont predykcji
D = N;      %horyzont dynamiki
Nu = 10;     %horyzont sterowania
lambda = 0.2; %kara za zmiany sterowania

M_p = zeros(N, D-1);
M = zeros(N, Nu);
delta_up = zeros(1, D-1);
t = zeros(1, k_end);
J = zeros(1, k_end);
e = 0;
delta_uk = 0;

%% Step response
u1 = zeros(1, N);
u2 = zeros(1, N);
u3 = zeros(1, N);
u4 = zeros(1, N);
u5 = zeros(1, N);
s1 = zeros(1, N);
s2 = zeros(1, N);
s3 = zeros(1, N);
s4 = zeros(1, N);
s5 = zeros(1, N);

%% Współczynniki wynikające z G_z
a2 = 0.02804;
a1 = 0.05607;
a0 = 0.02804;
b1 = 0.9592;
b0 = 0.1461;

u1(1:N) = -20;
u2(1:N) = -10;
u3(1:N) = 1;
u4(1:N) = 10;
u5(1:N) = 20;
s1(tau/Tp+1) = a2 * u1(1);
s1(tau/Tp+2) = b1*s1(tau/Tp+1) + a2 * u1(2) + a1 * u1(1);
s2(tau/Tp+1) = a2 * u2(1);
s2(tau/Tp+2) = b1*s2(tau/Tp+1) + a2 * u2(2) + a1 * u2(1);
s3(tau/Tp+1) = a2 * u3(1);
s3(tau/Tp+2) = b1*s3(tau/Tp+1) + a2 * u3(2) + a1 * u3(1);
s4(tau/Tp+1) = a2 * u4(1);
s4(tau/Tp+2) = b1*s4(tau/Tp+1) + a2 * u4(2) + a1 * u4(1);
s5(tau/Tp+1) = a2 * u5(1);
s5(tau/Tp+2) = b1*s5(tau/Tp+1) + a2 * u5(2) + a1 * u5(1);


for k=(tau/Tp+3):N
    s1(k) = b1 * s1(k-1) - b0 * s1(k-2) + a2 * u1(k-tau/Tp) + a1 * u1(k-1-tau/Tp) + a0 * u1(k-2-tau/Tp);
    s2(k) = b1 * s2(k-1) - b0 * s2(k-2) + a2 * u2(k-tau/Tp) + a1 * u2(k-1-tau/Tp) + a0 * u2(k-2-tau/Tp);
    s3(k) = b1 * s3(k-1) - b0 * s3(k-2) + a2 * u3(k-tau/Tp) + a1 * u3(k-1-tau/Tp) + a0 * u3(k-2-tau/Tp);
    s4(k) = b1 * s4(k-1) - b0 * s4(k-2) + a2 * u4(k-tau/Tp) + a1 * u4(k-1-tau/Tp) + a0 * u4(k-2-tau/Tp);
    s5(k) = b1 * s5(k-1) - b0 * s5(k-2) + a2 * u5(k-tau/Tp) + a1 * u5(k-1-tau/Tp) + a0 * u5(k-2-tau/Tp);
end