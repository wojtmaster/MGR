function [] = tf_discrete(G_z, u, t, k)

%% Równania dla dyskretnych równań stanu
h = lsim(G_z, u, t);

figure;
hold on;
stairs(0:k-1, h(:,1), 'r', 'LineWidth', 2);
hold off;
legend('h_1^{G(z)}');
xlabel('k');
ylabel('h [cm]');
title('Wysokość słupa cieczy w zbiorniku 1. - h_1(k)');
grid on;

figure;
hold on;
stairs(0:k-1, h(:,2), 'r', 'LineWidth', 2);
hold off;
legend('h_2^{G(z)}');
xlabel('k');
ylabel('h [cm]');
title('Wysokość słupa cieczy w zbiorniku 2. - h_2(k)');
grid on;
end