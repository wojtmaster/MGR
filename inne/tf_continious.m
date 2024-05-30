function [] = tf_continious(G_s, u, t)

%% Składowe odpowiedzi na poszczególne wymuszenia
h = lsim(G_s, u, t);

figure;
hold on;
plot(t, h(:,1), 'g', 'LineWidth', 2);
hold off;
legend('h_1^{G(s)}');
xlabel('t [s]');
ylabel('h [cm]');

figure;
hold on;
plot(t, h(:,2), 'g', 'LineWidth', 2);
hold off;
legend('h_2^{G(s)}');
xlabel('t [s]');
ylabel('h [cm]');
end