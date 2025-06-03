clear;
close;
clc;
%% Static characteristic
v = linspace(-2,2,100);
y = zeros(size(v));

for i = 1:length(v)
    y(i) = (0.3163*v(i)) / (sqrt(0.1 + 0.9*v(i)^2));
end

%% FIS
fis = sugfis("Name", "SugenoSystem");
fis = addInput(fis, [-2 2], "Name", "v_k");

% Dodanie dwch funkcji przynalenoci typu sigmoid
fis = addMF(fis, "v_k", "sigmf", [-5 0], "Name", "M1");  % malejc
fis = addMF(fis, "v_k", "sigmf", [5 0], "Name", "M2");   % rosnca

% plotmf(fis, 'input', 1);

% Dodanie wyjcia y_kplus1
fis = addOutput(fis, [-1 1], "Name", "y_k");

% Dodanie dwch funkcji wyjciowych typu constant
fis = addMF(fis, "y_k", "constant", -0.3289, "Name", "output_M1");
fis = addMF(fis, "y_k", "constant", 0.3289, "Name", "output_M2");

% Dodanie regu
rule1 = "if v_k is M1 then y_k is output_M1";
rule2 = "if v_k is M2 then y_k is output_M2";
fis = addRule(fis, [rule1; rule2]);

% Symulacja i wykres
v = linspace(-2, 2, 100)';
y_fuzzy = evalfis(fis, v);

% %% Wizualizacja
% figure;
% plot(v, y, 'b');
% hold on;
% plot(v, y_fuzzy, 'r')
% xlabel('v(k)');
% ylabel('y(k+1)');
% title('Wyjcie systemu Sugeno');
% legend('y', ' y_{fuzzy}');
% grid on;

%% DMC
N = 30;
D = 30;
Nu = 15;
lambda = 4;

%% Step respone
s = zeros(1, D);
u = ones(1, D);
a = [-1.4138 0.6065];
b = [0.1044 0.0883];
for k = 1:D
    if (k == 1)
        s(k) = b(1);
    elseif (k == 2)
        s(k) = -a(1)*s(k-1) + b*[u(k-1) u(k-1)]';
    else
        s(k) = -a*s(k-1:-1:k-2)' + b*u(k-1:-1:k-2)';    
    end
end
b(2)=0;
%% DMC - kontynuacja
M = zeros(N, Nu);
M_p = zeros(N, D-1);
            
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
ke = sum(K(1, :));          

% Implementacja macierzy M_p
for i = 1:N
    for j = 1:D-1
        if(i+j <= D)
            M_p(i,j) = s(i+j) - s(j);
        else 
            M_p(i,j) = s(D) - s(j);
        end
    end
end
ku = K(1, :) * M_p;

%% Sterowanie
clear v;
clear y;
clear u;

kk = 250;
y_zad = zeros(1,kk);
y_zad(1:kk/2) = 0.3;
y_zad(kk/2+1:end) = 0;

v.lmpc = zeros(1, kk);
v.fmpc = zeros(1, kk);
v.nmpc = zeros(1, kk);
u.lmpc = zeros(1, kk);
u.fmpc = zeros(1, kk);
u.nmpc = zeros(1, kk);
y.lmpc = zeros(1, kk);
y.fmpc = zeros(1, kk);
y.nmpc = zeros(1, kk);

delta_up.lmpc = zeros(1, D-1);
delta_up.fmpc = zeros(1, D-1);
delta_up.nmpc = zeros(1, D-1);
delta_uk.lmpc = 0;
delta_uk.fmpc = 0;
delta_uk.nmpc= 0;

dv = 0.001;

disp(b);

% Algorytm LMPC
for k = 3:kk
    v.lmpc(k) = -a*v.lmpc(k-1:-1:k-2)' + b*u.lmpc(k-1:-1:k-2)';

    y.lmpc(k) = 0.3163*v.lmpc(k) / sqrt(0.1 + 0.9*v.lmpc(k)^2);

    delta_up.lmpc = [delta_uk.lmpc, delta_up.lmpc(1:end-1)];

    e.lmpc = y_zad(k) - y.lmpc(k);

    delta_uk.lmpc = ke * e.lmpc - ku * delta_up.lmpc';

    u.lmpc(k) = u.lmpc(k-1) + delta_uk.lmpc;
end

% Algorytm FMPC
for k = 3:kk
    v.fmpc(k) = -a*v.fmpc(k-1:-1:k-2)' + b*u.fmpc(k-1:-1:k-2)';

    y.fmpc(k) = 0.3163*v.fmpc(k) / sqrt(0.1 + 0.9*v.fmpc(k)^2);
    % dydv = (y.fmpc(k) - evalfis(fis, (v.fmpc(k)-dv))) / dv;
    dydv = (evalfis(fis, (v.fmpc(k))) - evalfis(fis, (v.fmpc(k)-dv))) / dv;

    M_new = dydv*M;
    K_new = ((M_new' * M_new + lambda * eye(Nu))^(-1)) * M_new';

    y_0 = fuzzy_free_response(v.fmpc, u.fmpc, y.fmpc, k, a, b, N, fis);
    delta_uk.fmpc = K_new(1,:) * (y_zad(k) * ones(N,1) - y_0);
    u.fmpc(k) = u.fmpc(k-1) + delta_uk.fmpc;
end

% Algorytm NMPC
for k = 3:kk
    v.nmpc(k) = -a*v.nmpc(k-1:-1:k-2)' + b*u.nmpc(k-1:-1:k-2)';
    y.nmpc(k) = 0.3163*v.nmpc(k) / sqrt(0.1 + 0.9*v.nmpc(k)^2);
    cost_func = @(du) (sum((y_zad(k)*ones(N,1) - nonlinear_func(v.nmpc, u.nmpc, du, k, a, b, N)).^2) + lambda * sum(du.^2));
  
    % Ograniczenia na przyrosty sterowania
    lb = -0.1 * ones(Nu, 1);
    ub = 0.1 * ones(Nu, 1);
    
    % Optymalizacja
    options = optimoptions('fminunc', 'Display', 'off', ...
        'OptimalityTolerance', 1e-6, ...
        'StepTolerance', 1e-6);
    % delta_u = fmincon(cost_func, zeros(Nu, 1), [], [], [], [], [], [], [], options);
    delta_u = fminunc(cost_func, [delta_up.nmpc(1), delta_up.nmpc(1:end-1)], options);
    
    delta_up.nmpc = [delta_u(1), delta_up.nmpc(1:end-1)];
    
    % Aktualizacja sterowania
    u.nmpc(k) = u.nmpc(k-1) + delta_u(1);
end

figure;
subplot(1,2,1);
plot(0:kk-1, y_zad, 'k');
hold on;
plot(0:kk-1, y.lmpc, 'r');
plot(0:kk-1, y.fmpc, 'b');
plot(0:kk-1, y.nmpc, 'g');
xlabel('k');
ylabel('y(k)');
legend('y_{zad}', 'y', 'y_{wiener}', 'y_{nmpc}');
grid on;

subplot(1,2,2);
stairs(0:kk-1, u.lmpc, 'r');
hold on;
stairs(0:kk-1, u.fmpc, 'b.');
stairs(0:kk-1, u.nmpc, 'g-');
xlabel('k');
ylabel('u(k)');
legend('u', 'u_{wiener}', 'u_{nmpc}');
grid on;

function y_0 = fuzzy_free_response(v, u, y, k, a, b, N, fis)
    v_0 = zeros(N, 1);
    y_0 = zeros(N, 1);
    u_0 = u(k-1);

    for i = 1:N
        if i == 1
            v_0(i) = -a(1)*v(k) - a(2)*v(k-1) + b(1)*u(k-1) + b(2)*u(k-1);
        elseif i == 2
            v_0(i) = -a(1)*v_0(i-1) - a(2)*v(k) + b(1)*u_0 + b(2)*u(k-1);
        else
            v_0(i) = -a(1)*v_0(i-1) - a(2)*v_0(i-2) + b(1)*u_0 + b(2)*u_0;
        end
    end
    % v_0(1)=[];

    % Zakcenie DMC
    y_hat = evalfis(fis, v(k));  % biece przewidywane wyjcie
    dk = y(k) - y_hat;

    % Odpowied swobodna
    for i = 1:N
        y0 = evalfis(fis, v_0(i));
        y_0(i) = y0 + dk;             % rozmyta predykcja + zakcenie
    end
end


function y_pred = nonlinear_func(v, u, du, k, a, b, N)

    % Predykcja
    y_pred = zeros(N,1);
    u_pred = zeros(N,1);   % przysze wartoci sterowania
    v_pred = zeros(N,1);   % przysze wartoci stanu

    for i = 1:N
        % aktualizacja sterowania na horyzoncie sterowania, reszta jako
        % ostatni przyrost
        if i <= length(du)
            du_i = du(i);
        else
            du_i = du(end);
        end   

        % sterowanie w i-tym kroku predykcji
        if i == 1
            u_pred(i) = u(k-1) + du_i;
        else
            u_pred(i) = u_pred(i-1) + du_i;
        end

        % Symulacja modelu nieliniowego
        if i == 1
            v_pred(i) = -a(1)*v(k) - a(2)*v(k-1) + b(1)*u(k-1) + b(2)*u(k-1);
        elseif i == 2
            v_pred(i) = -a(1)*v_pred(i-1) - a(2)*v(k) + b(1)*u_pred(i-1) + b(2)*u(k-1);
        else
            v_pred(i) = -a(1)*v_pred(i-1) - a(2)*v_pred(i-2) + b(1)*u_pred(i-1) + b(2)*u_pred(i-2);
        end

        % if i>1
        %     % Nieliniowe wyjcie - predykcja
            y_pred(i) = (0.3163 * v_pred(i)) / sqrt(0.1 + 0.9 * v_pred(i)^2);
        % end
    end
end