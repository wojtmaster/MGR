function [K] = GPC_matrix(N, Nu, lambda, s)

    M = zeros(N, Nu);
    
    for i = 1:N
        for j = 1:Nu
            if(i >= j)
                M(i, j) = s(i-j+1);
            end
        end
    end
    
    K = ((M' * M + lambda * eye(Nu))^(-1)) * M';

end