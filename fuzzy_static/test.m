%% DANE
clear;

a = 540;
c = 0.85;
alpha_1 = 26;
alpha_2 = 20;

% Punkt pracy
F_10 = 90;
F_D0 = 30;
h_10 = ((F_10+F_D0)/alpha_1)^2;
h_20 = ((F_10+F_D0)/alpha_2)^2;
V_10 = a * h_10;
V_20 = c * h_20^2;

% Opóźnienie
tau = 100;
% Okres próbkowania
Tp = 20;
% Próbki dyskretne
kk = 500;
% Wektor czasu
t = 0:Tp:(kk-1)*Tp;

[G_s, G_z] = tf_function(a, c, alpha_1, alpha_2, V_10, V_20, tau, Tp);

b(1) = [G_z.Numerator{2,1}(2)];
b(2) = [G_z.Numerator{2,1}(3)];

a(1) = [G_z.Denominator{2,1}(2)];
a(2) = [G_z.Denominator{2,1}(3)];

u = linspace(0, 180, 1000);
y = zeros(size(u));

for i = 1:length(u)
    y(i) = ((u(i)+F_D0)/alpha_2)^2;
end

R = cell(1,2);
alpha_1 = 1;

% figure; 
% plot(u, y);