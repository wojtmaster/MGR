function [M, M_p, K, ke, ku] = macierze_DMC(N, D, Nu, lambda, s)

%% Skrypt generujący macierze M, M_p, K dla zadanych wartości N, Nu, D, lambda oraz rzędnych odp. skokowej

M = zeros(N, Nu);
M_p = zeros(N, D-1);

%% Implementacja macierzy M_p
for i = 1:N
    for j = 1:D-1
        if(i+j <= D)
            M_p(i, j) = s(i+j) - s(j);
        else 
            M_p(i, j) = s(D) - s(j);
        end
    end
end

%% Implementacja macierzy M
for i = 1:N
    for j = 1:Nu
        if(i >= j)
            M(i, j) = s(i-j+1);
        end
    end
end

%% Macierz K
K = ((M' * M + lambda * eye(Nu))^(-1)) * M';
ke = sum(K(1, :));
ku = K(1, :) * M_p;

end