%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kaskadowa regulacja poziomu wody w zbiornikach %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Rozmyta regulacja predykcyjna DMC na podstawie odpowiedzi skokowej na różne wymuszenia skokowe
% Autor: Wojciech Rogalski
% Data: 25.10.2023r.
clear all;

%% DMC parameters
N = 100;         %horyzont predykcji
D = N;           %horyzont dynamiki
Nu = 5;          %horyzont sterowania
lambda = 1;      %kar za zmiany sterowania

%% Dane
steps = 500;    %liczba kroków symulacji
tau = 100;      %bezwładność układu
Tp = 50;        %okres próbkowania

t = zeros(1, steps);
error = zeros(1, steps);

%% Romzywanie u
sets = 5;
width = 10;
n = 1000;
u_Fmin = -90;
u_Fmax = 90;
u_F = zeros(1, sets);
u_F(1) = u_Fmin;
u_F(2) = -60;
u_F(3) = 1;
u_F(4) = 60;
u_F(5) = u_Fmax;

x = linspace(u_Fmin, u_Fmax, n);
w_arr = cell(1, sets);
w = zeros(1, sets);

for i = 1:sets
    if i == 3
        w_arr{i} = gaussmf(x, [1.5*width, u_F(i)]);
    else
        w_arr{i} = gaussmf(x, [width, u_F(i)]);
    end
end

figure(1);
hold on;
plot(x, w_arr{1});
plot(x, w_arr{2});
plot(x, w_arr{3});
plot(x, w_arr{4});
plot(x, w_arr{5});
hold off;
title('Funkcje przynależności wartości u');
xlabel('u');
ylabel('mi(u)');
axis([u_Fmin u_Fmax 0 1]);

%% Rzędne odpowiedzi skokowej dla poszczególnych wymuszeń
s_F = cell(1, sets);
M_F = cell(1, sets);
M_pF = cell(1, sets);
K_F = cell(1, sets);
ke_F = zeros(1, sets);
ku_F = cell(1, sets);

for i = 1:sets
    [s_F{i}, a2, a1, a0, b1, b0] = step_response_fuzzy(u_F(i), N, Tp); 
    [M_F{i}, M_pF{i}, K_F{i}, ke_F(i), ku_F{i}] = macierze_DMC(N, D, Nu, lambda, s_F{i});
end

%% Alokacja pamięci
u = zeros(1, steps);
y = zeros(1, steps);
y_zad = zeros(1, steps);

%% Zadane wartości sygnału regulowanego
y_zad(1:100) = 10;
y_zad(101:200) = -5;
y_zad(201:300) = 0;
y_zad(301:400) = 5;
y_zad(401:steps) = -10;

delta_ukF = zeros(1, sets);
delta_uk = 0;
delta_up = zeros(1, D-1);
e = 0;

%% Ograniczenia sygnału sterującego i wyjścia obiektu
u_max = 20;
delta_uk_max = 2;
y_max = 20;

%% Sterowanie DMC
for i = 2:steps
    % Mapowanie wartości u(i-1) na wartości zbióru rozmytego
    number = round((u(i-1) - u_Fmin) * (n - 0)/(u_Fmax - u_Fmin) + 0);
    % Ograniczenia na wypadek wyjścia poza przedział
    if number > 1000
        number = 1000;
    elseif number < 1
        number = 1;
    end
    
    % Dobór wag ze względu na otrzymany numer
    for j = 1:sets
        w(j) = w_arr{j}(number);
    end

    if(i < tau/Tp+1)
        y(i) = 0;
    elseif(i == tau/Tp+1)
        y(i) = a2 * u(i-tau/Tp);
    elseif(i == tau/Tp+2)
        y(i) = b1*y(i-1) + a2 * u(i-tau/Tp) + a1*u(i-1-tau/Tp);
    else
        y(i) = b1*y(i-1) - b0*y(i-2) + a2 * u(i-tau/Tp) + a1 * u(i-1-tau/Tp) + a0 * u(i-2-tau/Tp);
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

    delta_uk_sum = 0;
    w_sum = 0;
    
    % Dyfuzyfikacja
    for j = 1:sets
        delta_ukF(j) = ke_F(j) * e - ku_F{j} * delta_up';
        delta_uk_sum = delta_uk_sum + delta_ukF(j) * w(j);
        w_sum = w_sum + w(j);
    end
    
    % Obliczenie przyrostu sterowania dla chwili (i+1) w chwili i-tej
    delta_uk = delta_uk_sum / w_sum;
    
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
    error(i) = error(i-1) + 1/2 * (y_zad(i) - y(i))^2;
end

disp(sum(error));

%% Prezentacja wyników
figure(2);
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
axis([0 steps -20 20]);
hold off;
legend('u');
title('u(k)');
xlabel('k');
ylabel('u');