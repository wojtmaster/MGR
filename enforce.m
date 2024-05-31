function [u] = enforce(kk)
    u = zeros(2,kk);
    % u(1, :) = F_1
    u(1, 1 : 0.25*kk) = 0;
    u(1, 0.25*kk+1 : 0.5*kk) = 2;
    u(1, 0.5*kk+1 : 0.75*kk) = -2;
    u(1, 0.75*kk+1 : end) = 4;
    for i = 1:kk
        if(u(1,i) < -90)
            u(1,i) = -90;
        end
    end
    
    % u(2, :) = F_D
    u(2, 1 : 0.2*kk) = 2;
    u(2, 0.2*kk+1 : 0.4*kk) = 2;
    u(2, 0.4*kk+1 : 0.6*kk) = -2;
    u(2, 0.6*kk+1 : 0.8*kk) = 4;
    u(2, 0.8*kk+1 : end) = 0;
    for i = 1:kk
        if(u(2,i) < -30)
            u(2,i) = -30;
        end
    end
end