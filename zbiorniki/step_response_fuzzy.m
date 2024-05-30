function [s, a2, a1, a0, b1, b0] = step_response_fuzzy(U, N, Tp)
    
    %% Skrypt generujący odpowiedzi skokowe dla wymuszenia skokowego
    
    %% Dane
    P = 540;
    C = 0.85;
    alfa_1 = 26;
    alfa_2 = 20; 
    
    %% Punkt pracy
    F_10 = 110;
    F_D0 = 30;
    h_10 = ((F_10+F_D0)/alfa_1)^2;
    h_20 = ((F_10+F_D0)/alfa_2)^2;
    tau = 100;
    
    %% Równania stanu
    A = [-alfa_1/(2*P*sqrt(h_10)) 0; alfa_1/(4*C*h_20*sqrt(h_10)) -alfa_2/(4*C*h_20*sqrt(h_20))];
    B = [1/P 1/P; 0 0];
    C = [0 1];
    
    s1 = tf('s');
    G_s = C*(s1*eye(2) - A)^(-1)*B;
    G_z = c2d(G_s, Tp, 'tustin');
    
    %% Step response
    u = zeros(1, N);
    s = zeros(1, N);
    
    %% Współczynniki wynikające z G_z
    a2 = 0.02804;
    a1 = 0.05607;
    a0 = 0.02804;
    b1 = 0.9592;
    b0 = 0.1461;
    
    u(1:N) = U;
    s(tau/Tp+1) = a2 * u(1);
    s(tau/Tp+2) = b1*s(tau/Tp+1) + a2 * u(2) + a1 * u(1);
    
    for i=(tau/Tp+3):N
        s(i) = b1 * s(i-1) - b0 * s(i-2) + a2 * u(i-tau/Tp) + a1 * u(i-1-tau/Tp) + a0 * u(i-2-tau/Tp);
    end
end