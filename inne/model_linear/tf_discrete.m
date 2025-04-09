function [y] = tf_discrete(G_z, u, t, kk)

%% Równania dla dyskretnych równań stanu
y = lsim(G_z, u, t);

% figure;
% hold on;
% stairs(0:k-1, h(:,1), 'r', 'LineWidth', 2);
% hold off;
% legend('h_1^{G(z)}');
% xlabel('k');
% ylabel('h [cm]');
% title('Wysokość słupa cieczy w zbiorniku 1. - h_1(k)');
% grid on;

figure;
hold on;
stairs(0:kk-1, y, 'r', 'LineWidth', 2);
hold off;
xlabel('k');
ylabel('y(k)');
title('Wysokość słupa cieczy w zbiorniku 2. - y(k)');
grid on;
end