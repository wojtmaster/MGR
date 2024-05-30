    function [G_z] = tf_function(a, c, alpha_1, alpha_2, V_10, V_20, tau, Tp)

%% Funckja do przejścia od równań stanu do transmitancji
% [x1, x2] = [V_1, V_2]
% [y1, y2] = [h_1, h_2]
% [u1, u2, u3] = [F_1, F_D]

A = [-alpha_1/(2*sqrt(a*V_10)), 0;
     alpha_1/(2*sqrt(a*V_10)), -alpha_2/(4*c^(1/4)*V_20^(3/4))];
B = [1 1; 0 0];
C = [1/a 0; 0 1/(2*sqrt(V_20/c))];
D = [0 0; 0 0];

% sys = ss(A,B,C_matrix,D);
% G_s = tf(sys);
%% Dodanie opóźnienia
s = tf('s');
G_s = C*(s*eye(2) - A)^(-1)*B + D;
G_s(1,1) = tf(G_s.Numerator(1,1), G_s.Denominator(1,1), 'InputDelay', tau);
G_s(2,1) = tf(G_s.Numerator(2,1), G_s.Denominator(2,1), 'InputDelay', tau);

%% Transmitancja dyskretna
G_z = c2d(G_s, Tp, 'zoh');
G_z.Variable = 'z^-1';

%% Odpowiedź skokowa
% figure;
% step(G_s, G_z);
end