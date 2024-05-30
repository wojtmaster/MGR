clear;
% s = [0 0 0.2 0.5 0.6 0.62];
N = 6;
D = N;
Nu = 3;
lambda = 0.01;

na = 1;
nb = 4;

kk = 15;

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

for i = 1:N
    for j = 1:Nu
        if(i >= j)
            M(i, j) = s(i-j+1);
        end
    end
end

K = ((M' * M + lambda * eye(Nu))^(-1)) * M';
K_1 = K(1, :);



y_0 = zeros(N,1);
y = zeros(kk,1);
y_zad = zeros(kk+N,1);
y_zad(2:end) = 1;
Y_zad = zeros(N,1);
Y_0 = zeros(N,1);

u = zeros(1,N);
d = zeros(1, kk);

for i = 2:kk
    if i == 2 || i == 3
        y(i) = -a(1)*y(i-1);
    elseif i == 4
        y(i) = -a(1)*y(i-1) + b(3)*u(i-3);
    else
        y(i) = -a(1)*y(i-1) + b(3)*u(i-3) + b(4) * u(i-4);
    end
    
    for k = 1:N
        if i == 2 && k == 1
            y_0(k) = -a(1)*y(i);
        elseif i == 2 && k == 2
            y_0(k) = -a(1)*y_0(k-1) + b(3)*u(i-1);
        elseif i == 3 && k == 1
            y_0(k) = -a(1)*y(i) + b(3)*u(i-2);
        else    
            if k == 1
                y_0(k) = -a(1)*y(i) + b(3)*u(i-2) + b(4) * u(i-3);
            elseif k == 2
                y_0(k) = -a(1)*y_0(k-1) + b(3)*u(i-1) + b(4) * u(i-2);
            else
                y_0(k) = -a(1)*y_0(k-1) + b(3)*u(i-1) + b(4) * u(i-1);
            end
        end
    end

    Y_zad = y_zad(i:i+(N-1));
    Y_0 = y_0(1:N);
    delta_u = K_1*(Y_zad - Y_0);
    if delta_u > 0.5
        delta_u = 0.5;
    end
    u(i) = u(i-1) + delta_u;
    if u(i) > 2
        u(i) = 2;
    end
end

figure;
plot(1:kk, y);

figure;
stairs(1:kk, u);