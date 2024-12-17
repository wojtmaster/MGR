%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   SZAU - PROJEKT - ZADANIE P14    %%
%%   Hubert Pisiecki 303170          %%
%%   Bartłomiej Baścik 336341        %%
%%   SEMESTR ZIMOWY 2024/2025        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Algorytm regulacji predykcyjnej w wersji analitycznej z pomiarem zakłóceń - DMC

clear all;
clc;
close all;

dmc_fig = figure('Color', 'w');
%% Parametry DMC
D = 70;         %horyzont dynamiki
Dz = 70;
N = 70;         %horyzont predykcji
Nu = 14;        %horyzont sterowania
lambda = 0.5;   %kara

%% Dane
A_2 = 310;
C_1 = 0.65;
alpha_1 = 17;
alpha_2 = 22; 

steps = 500;
tau = 100;
Tp = 10;

t = zeros(1, steps);
error = zeros(1, steps);

%% Punkt pracy
F_10 = 75;
F_D0 = 15; 
h_10 = ((F_10+F_D0)/alpha_1)^2;
h_20 = ((F_10+F_D0)/alpha_2)^2;

%% Równania stanu
A = [-alpha_1/(4*C_1*h_10*sqrt(h_10)) 0;...
      alpha_1/(2*A_2*sqrt(h_10)) -alpha_2/(2*A_2*sqrt(h_20))];
B = [1/(2*h_10*C_1) 1/(2*h_10*C_1); 0 0];
C_1 = [0 1];
D_1 = [0 0];
    
s = tf('s');
G_s = C_1*(s*eye(2) - A)^(-1)*B;
G_z = c2d(G_s, Tp, 'tustin');

model__ss = ss(A,B,C_1,D_1,'InputDelay',[100 0]);
G_s_test = tf(model__ss);
G_s_test_z = c2d(G_s_test, Tp, 'tustin');
G_zaklocenie_tf = G_s_test(1,2);

%% Badanie sterowalności

if rank([B A*B]) == size(A, 1)
    disp('Układ jest sterowalny')
else
    disp('Układ nie jest sterowalny')
end

%% Współczynniki transmitancji a oaz b
a1 = G_z.Denominator{1}(2);
a0 = G_z.Denominator{1}(3);
b2 = G_z.Numerator{1}(1);
b1 = G_z.Numerator{1}(2);
b0 = G_z.Numerator{1}(3);

u = zeros(1, D);
s = zeros(1, D);

u(1:D) = 1;
s(tau/Tp+1) = b2 * u(1);
s(tau/Tp+2) = -a1*s(tau/Tp+1) + b2 * u(2) + b1 * u(1);
    
for i=(tau/Tp+3):D
s(i) = -a1 * s(i-1) - a0 * s(i-2) +...
        b2 * u(i-tau/Tp) + b1 * u(i-1-tau/Tp) + b0 * u(i-2-tau/Tp);
end

%% odp skokowa zaklocenia
timefinal = Dz * Tp;
S_z = DMC_zaklocenia_odp_skok_Sz(Tp,G_zaklocenie_tf,timefinal);

%% K
M_p = oblicz_Mp(N,D,s);
M = oblicz_M(N,Nu,s);

%zaklocenia:
M_zp = DMCmatrixMzP(S_z,N,Nu);

K = ((M' * M + lambda * eye(Nu))^(-1)) * M';
K1 = K(1,:);
ke = sum(K1);
ku = K1 * M_p;
kz = K1 * M_zp;

u = zeros(1, steps);
y_Fin = zeros(1,steps);
y_Fd = zeros(1,steps);
y = zeros(1, steps);

%% Zadana wartość wielości regulowanej h_2
h_2zad = zeros(1, steps);
h_2zad(1:100) = 0;
h_2zad(101:300) = 50;
h_2zad(301:450) = 25;
h_2zad(450:steps) = 0;
%h_2zad(101:200) = 10;
%h_2zad(201:300) = 5;
%h_2zad(301:400) = -5;
%h_2zad(401:steps) = 0;

%% Trajektoria zaklocen:
Fd = zeros(1,steps);
Fd(1:50) = F_D0 ;
Fd(51:200) = 0;
Fd(201:300) = 2*F_D0;
Fd(301:500) = F_D0;

delta_up = zeros(1, D-1);
delta_z = zeros(1, D);
delta_uk = 0;
e = 0;

%% Ograniczenia sygnału sterującego i sterowanego
u_max = 500; %bylo 50
u_min = -100; %bylo 0
delta_uk_max = 150; %bylo 15
y_max = 300; %bylo 30
y_min = 0;
%% Algorytm DMC
for i = 2:steps
    if(i < tau/Tp+1)
        y_Fin(i) = 0;
    elseif(i == tau/Tp+1)
        y_Fin(i) = b2 * u(i-tau/Tp);
    elseif(i == tau/Tp+2)
        y_Fin(i) = -a1*y_Fin(i-1) + b2*u(i-tau/Tp) + b1*u(i-1-tau/Tp);
    else
        y_Fin(i) = -a1*y_Fin(i-1) - a0*y_Fin(i-2) +...
                    b2*u(i-tau/Tp) + b1*u(i-1-tau/Tp) + b0*u(i-2-tau/Tp);
    end

    if i < 1
        y_Fd(i)
    elseif i == 1
        y_Fd(i) = b2*Fd; 
    elseif i == 2
        y_Fd(i) = -a1*y_Fd(i-1) + b2*Fd(i) + b1*Fd(i-1); 
    else
        y_Fd(i) = -a1*y_Fd(i-1) - a0*y_Fd(i-2) +...
                   b2*Fd(i) + b1*Fd(i-1) + b0*Fd(i-2); 
    end 

    y(i) = y_Fin(i) + y_Fd(i);

    if(y(i) > y_max)
        y(i) = y_max;
    elseif(y(i) < y_min)
        y(i) = -y_min;
    end
   
        for j = (D-1):-1:1
            if j == 1
                delta_up(j) = delta_uk;
            else
                delta_up(j) = delta_up(j-1);
            end
        end 
    %Obliczenie zmiany zaklocenia:
    delta_z = [Fd(i) - Fd(i-1) , delta_z(1,1:(end-1))];


    % Błąd sterowania 
    e(i) = h_2zad(i) - y(i);
    
    % Obliczenie przyrostu sterowania dla chwili (i+1) w chwili i-tej
    % Z POMIAREM ZAKLOCEN:
    delta_uk = ke * e(i) - ku * delta_up' - kz * delta_z';
    % BEZ POMIARU ZAKLOCEN:
    % delta_uk = ke * e(i) - ku * delta_up';
    
    if(delta_uk > delta_uk_max)
        delta_uk = delta_uk_max;
    elseif(delta_uk < -delta_uk_max)
        delta_uk = -delta_uk_max;
    end
    
    % Obliczenie sterowania dla chwili (i+1) w chwili i-tej
    u(i) = u(i-1) + delta_uk; 
    
    if(u(i) > u_max)
        u(i) = u_max;
    elseif(u(i) < u_min)
        u(i) = u_min;
    end

    t(i) = (i-1)*Tp;

    % Obliczenie wskaźnika regulacji
    error(i) = error(i-1) + 1/2 * (h_2zad(i) - y(i))^2;
end
disp(sum(error));

%% Generowanie wykresów
figure(dmc_fig);
subplot(2,1,1);
grid on;
hold on;
stairs(y);
stairs(h_2zad);
hold off;
title('h_{2}(k)');
xlabel('k');
ylabel('h_{2}');
% legend('bez pomiaru zakl.','wartosc zadana', 'z pomiarem zakl', 'Location','northeast');
% axis([0 steps -10 20]);

subplot(2,1,2);
grid on;
hold on;
stairs(u);
% axis([0 steps -5 50]);
hold off;
legend('F_{1in}');
title('F_{1in}(k)');
xlabel('k');
ylabel('F_{1in}');
% legend('bez pomiaru zakl.', 'z pomiarem zakl', 'Location','northeast');
% exportgraphics(gcf,'..\..\figures\zad1\DMC\DMC_lambda_5.pdf','Resolution',400)