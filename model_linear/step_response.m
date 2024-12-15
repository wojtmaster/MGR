function [s] = step_response(A, B, na, nb, N)
    s = zeros(1,N);
    b = zeros(1,nb);
    a = A;

    b(6) = B(1);
    b(7) = B(2);

    s(1) = b(1);
    s(2) = sum(b(1:2)) - a(1)*s(1); 
    for k = 3:N
        s(k) = sum(b(1:min(k,nb))) - a(1:min(k-1, na))*[s(k-1) s(k-na)]';
    end
end