function [a, b, K] = sopdt(y, t, tau, Tp)
    T_0 = tau;
    K_0 = 0.6025;
    T_1 = 212;
    T_2 = 15;
    num = K_0;
    den =conv([T_1 1], [T_2 1]);
    
    G_s = tf(num, den, 'InputDelay', T_0);
    
    G_z = c2d(G_s, Tp, 'zoh');
    G_z.Variable = 'z^-1';
    
    a(1:2) = G_z.Denominator{1}(2:end);
    b(1:2) = G_z.Numerator{1}(2:end);
    K = dcgain(G_z);

    %% Prezentacja wyników SOPDT
    figure;
    plot(t, y, 'r-','LineWidth',2);
    hold on;
    step(G_z);
    grid on;
    legend('h', 'h^{mod}', 'Location', 'northwest');
    % file_name = sprintf('../raport/pictures/model_sopdt.pdf');
    % exportgraphics (gcf, file_name);
end