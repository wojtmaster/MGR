%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kaskadowa regulacja poziomu wody w zbiornikach %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Regulacja predykcyjna DMC na podstawie odpowiedzi skokowej na skokowe wymuszenie jednostkowe
% Autor: Wojciech Rogalski
% Data: 25.10.2023r.
clear all;

%% DMC parameters
N = 100;        %horyzont predykcji
D = N;          %horyzont dynamiki
Nu = 5;         %horyzont sterowania
lambda = 1;     %kara za zmiany sterowania

%% Dane
steps = 500;    %liczba kroków symulacji
tau = 100;      %bezwładność układu
Tp = 50;        %okres próbkowania

t = zeros(1, steps);
error = zeros(1, steps);

%% Odpowiedzi na wymuszenie jednostkowe
[s, a2, a1, a0, b1, b0] = step_response(N, Tp);

%% Implementacja macierzy M_p, M, K
[M, M_p, K, ke, ku] = macierze_DMC(N, D, Nu, lambda, s);

%% Alokacja pamięci
u = zeros(1, steps);
y = zeros(1, steps);
y_zad = zeros(1, steps);

%% Zadane wartości sygnału regulowanego
y_zad(1:100) = 20;
y_zad(101:200) = 10;
y_zad(201:300) = 0;
y_zad(301:400) = -10;
y_zad(401:steps) = -20;

delta_up = zeros(1, D-1);
delta_uk = 0;
e = 0;

%% Ograniczenia sygnału sterującego i wyjścia obiektu
u_max = 20;
delta_uk_max = 5;
y_max = 20;

%% Sterowanie DMC
for i = 2:steps
    if(i < tau/Tp+1)
        y(i) = 0;
    elseif(i == tau/Tp+1)
        y(i) = a2 * u(i-tau/Tp);
    elseif(i == tau/Tp+2)
        y(i) = b1*y(i-1) + a2 * u(i-tau/Tp) + a1*u(i-1-tau/Tp);
    else
        y(i) = b1*y(i-1) - b0*y(i-2) + a2 * u(i-tau/Tp) + a1*u(i-1-tau/Tp) + a0 * u(i-2-tau/Tp);
    end
    
    % Ograniczenia wartości sygnału wyjściowego, tj. wysokości h_2
    if(y(i) > y_max)
        y(i) = y_max;
    elseif(y(i) < -y_max)
        y(i) = -y_max;
    end

    % Przepisanie sterowań do wektora przeszłych sterowań
        for j = (D-1):-1:1
            if j == 1
                delta_up(j) = delta_uk;
            else
                delta_up(j) = delta_up(j-1);
            end
        end 
    % Oblicznie uchybu    
    e = y_zad(i) - y(i);

    % Obliczenie przyrostu sterowania dla chwili (i+1) w chwili i-tej
    delta_uk = ke * e - ku * delta_up';
    
    % Ograniczenie wartości przyrostu sterowania
    if(delta_uk > delta_uk_max)
        delta_uk = delta_uk_max;
    elseif(delta_uk < -delta_uk_max)
        delta_uk = -delta_uk_max;
    end
    
    % Obliczenie sterowania dla chwili (i+1) w chwili i-tej
    u(i) = u(i-1) + delta_uk; 
    
    % Ograniczenie sterowania
    if(u(i) > u_max)
        u(i) = u_max;
    elseif(u(i) < -u_max)
        u(i) = -u_max;
    end

    t(i) = (i-1)*Tp;

    % Obliczenie wskaźnika regulacji
    error(i) = error(i-1) + 1/2 * (y_zad(i) - y(i))^2;
end
disp(sum(error));

%% Prezentacja wyników
subplot(2,1,1);
hold on;
stairs(y);
stairs(y_zad);
hold off;
legend('h_{2}', 'h_{2zad}');
title('h_{2}(k)');
xlabel('k');
ylabel('h_{2}');
axis([0 steps -20 20]);

subplot(2,1,2);
hold on;
stairs(u);
axis([0 steps -40 40]);
hold off;
legend('u');
title('u(k)');
xlabel('k');
ylabel('u');