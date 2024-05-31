% function [] = GPC_analitic()
clear;
% s = [0 0 0.2 0.5 0.6 0.62];
N = 6;
D = N;
Nu = 3;
lambda = 5;

Tp = 1;
% Próbki dyskretne
kk = 500;
% Wektor czasu
t = 0:Tp:(kk-1)*Tp;

na = 2;
nb = 4;

a(1) = -0.2676;
b(1:2) = 0;
b(3) = 0.1989;
b(4) = 0.2552;

s(1) = b(1);
s(2) = sum(b(1:2)) - a(1)*s(1); 
for k = 3:N
    s(k) = sum(b(1:min(k,nb))) - a(1)*s(k-1);
end

M = zeros(N, Nu);

for k = 1:N
    for i = 1:Nu
        if(k >= i)
            M(k, i) = s(k-i+1);
        end
    end
end

K = ((M' * M + lambda * eye(Nu))^(-1)) * M';
K_1 = K(1, :);

% Dane równania
A = [1 -1.2676 0.2676]; % współczynniki dla A(z^-1)
B = [0.1989 0.2552]; % współczynniki dla B(z^-1)

tau = 2;

E = cell(N,1);
F = cell(N,1);
G = cell(N,1);
G_PG = cell(N-tau, 1);

for k = 1:N
    E{k} = zeros(1,k);
    F{k} = zeros(1,na);
    G{k} = zeros(1, k+1);
end

E{1} = 1;
F{1} = -A(2:end);

for p = 2:N
    % E
    for k = 1:p
        if k <= p-1
            e = E{p-1}(k);
        else
            e = F{p-1}(1);
        end
        E{p}(k) = e;
    end

    % F
    for k = 1:na
        if k < na
            f = -F{p-1}(1) * A(k+1) + F{p-1}(k+1);
        else
            f = -F{p-1}(1) * A(k+1);
        end
        F{p}(k) = f;
    end
end
% opóźnienie tau = 2 --> zaczynamy od p = 3
for p = tau+1:N
    for k = 1:length(E{p}) + 1
        if k == 1
            G{p}(k) = E{p}(1)*B(1);
        elseif k < length(E{p}) + 1
            G{p}(k) = E{p}(k-1)*B(2) + E{p}(k)*B(1);
        else
            G{p}(k) = E{p}(k-1)*B(2);
        end
    end
end

for k = 1:N-tau
    G_PG{k} = G{k+tau}(end-(Nu-1):end);
end
G_PG = cell2mat(G_PG);
F = cell2mat(F);
k_y = zeros(1,na);
k_u = zeros(1,nb-1);
for k = 1:na
    k_y(k) = K_1*F(:,k);
end

for k = 1:nb-1
    k_u(k) = K_1(tau+1:end)*G_PG(:,k);
end
k_e = sum(K_1);

y = zeros(kk,1);
y_zad = zeros(kk,1);
y_zad(2:end) = 1;

delta_up = zeros(1, nb-1);
delta_u = 0;

u = zeros(1,kk);
u_max = 2;
delta_uk_max = 0.5;

for k = 2:kk
    if k == 2 || k == 3
        y(k) = -a(1)*y(k-1);
    elseif k == 4
        y(k) = -a(1)*y(k-1) + b(3)*u(k-3);
    else
        y(k) = -a(1)*y(k-1) + b(3)*u(k-3) + b(4) * u(k-4);
    end

    for i = (nb-1):-1:1
        if i == 1
            delta_up(i) = delta_u;
        else
            delta_up(i) = delta_up(i-1);
        end
    end
    delta_up = [delta_u, delta_up(1:end-1)];

    if k == 2
        delta_u = k_e*y_zad(k) - k_y*[y(k-1); 0] - k_u*delta_up';
    else
        delta_u = k_e*y_zad(k) - k_y*[y(k-1); y(k-2)] - k_u*delta_up';
    end
    if delta_u > delta_uk_max
        delta_u = delta_uk_max;
    elseif delta_u < -delta_uk_max
        delta_u = -delta_uk_max;
    end
    u(k) = u(k-1) + delta_u;
    if u(k) > u_max
        u(k) = u_max;
    elseif u(k) < -u_max
        u(k) = -u_max;
    end
end

figure;
hold on;
plot(1:kk, y);
stairs(1:kk, y_zad);
hold off;

figure;
stairs(0:(kk-1), u);
% end

%%
% GPC ale z delta_u = k_e*y_zad(k) - k_y*[y(k-1); y(k-2)] - k_u*delta_up';

function [] = test(a, b, N, D, Nu, lambda, y_zad, s, kk, t, tau)
    M = zeros(N, Nu);
    
    % Implementacja macierzy M
    for i = 1:N
        for j = 1:Nu
            if(i >= j)
                M(i,j) = s(i-j+1);
            end
        end
    end

    % Macierz K
    K = ((M' * M + lambda * eye(Nu))^(-1)) * M';
    K_1 = K(1,:);
    
    % Alokacja pamięci
    y = zeros(1, kk);
    delta_up = zeros(D-1,1);
    delta_u = 0;
    u = zeros(1, kk);
    
    % Sterowanie DMC
    for k = 2:kk
        if k < tau+2
            y(k) = 0;
        elseif k == tau+2
            y(k) = - a(1)*y(k-1) - a(2)*y(k-2) + b(1)*u(k-(tau+1));
        else
            y(k) = - a(1)*y(k-1) - a(2)*y(k-2) + b(1)*u(k-(tau+1)) + b(2)*u(k-(tau+2));
        end
        
        % Przepisanie sterowań do wektora przeszłych sterowań
            for i = (D-1):-1:1
                if i == 1
                    delta_up(i) = delta_u;
                else
                    delta_up(i) = delta_up(i-1);
                end
            end 

        Y = ones(N,1)*y(k);
        Y_0 = Y + M_p*delta_up;
        Y_zad = ones(N,1)*y_zad(k);

        % Obliczenie przyrostu sterowania dla chwili (i+1) w chwili i-tej
        delta_u = K_1*(Y_zad - Y_0);
        
        % Obliczenie sterowania dla chwili (i+1) w chwili i-tej
        u(1,k) = u(1,k-1) + delta_u; 
    end
    
    % Prezentacja wyników
    figure;
    hold on;
    stairs(1:kk, y_zad, 'r');
    stairs(1:kk, y, 'b');
    hold off;
    legend('y_{zad}', 'y');
    plot_title = sprintf('Sygnał wyjściowy y(k) \n Regulator DMC');
    title(plot_title);
    xlabel('k');
    ylabel('y(k)');
    
    figure;
    hold on;
    stairs(0:kk-1, u, 'r');
    hold off;
    legend('u');
    plot_title = sprintf('Sygnał sterujący u(k) \n Regulator DMC');
    title(plot_title);
    xlabel('k');
    ylabel('u');
end