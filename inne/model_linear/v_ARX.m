function [v] = v_ARX(a, b, K, u, y, kk, tau)

    v = zeros(1, kk);
    v(1:tau+2) = y(1:tau+2);

    % Model ARX
    % Zbiór uczący
    for k = tau+3:kk
        v(k) = - a'*[y(k-1:-1:k-2)]' + (b(tau+1:tau+2)/K)'*[u(1,k-(tau+1):-1:k-(tau+2))]' + b(tau+1:tau+2)'*[u(2, k-1:-1:k-2)]';
    end 

end