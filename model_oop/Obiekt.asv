classdef Obiekt
    properties
        % Stałe
        A = 540;
        C = 0.85;
        alpha_1 = 26;
        alpha_2 = 20;
        
        % Punkt pracy
        F_10 = 90;
        F_D0 = 30;
        
        % Opóźnienie
        tau = 100;
        % Okres próbkowania
        Tp = 20;
        % Próbki dyskretne
        kk = 5000;
        % Podział danych
        n = 500;
        
        % Do obliczenia
        h_10
        h_20
        V_10
        V_20

        % Mode - ARX / OE
        mode
        % Następniki - linear / nonlinear
        s
        % Współczynniki modelu
        delay
    end

    methods
        function obj = Obiekt(mode, s)
            obj.h_10 = ((obj.F_10+obj.F_D0)/obj.alpha_1)^2;
            obj.h_20 = ((obj.F_10+obj.F_D0)/obj.alpha_2)^2;
            obj.V_10 = obj.A * obj.h_10;
            obj.V_20 = obj.C * obj.h_20^2;

            obj.mode = mode;
            obj.s = s;
            obj.delay = obj.tau/obj.Tp;
        end

        function [a, b, K, G_z] = sopdt(obj, K_0, T_0, T_1, T_2)
            %% Model inercjalny SOPDT na podstawie odpwiedzi skokowej
            num = K_0;
            den = conv([T_1 1], [T_2 1]);
            
            G_s = tf(num, den, 'InputDelay', T_0);
            
            G_z = c2d(G_s, obj.Tp, 'zoh');
            G_z.Variable = 'z^-1';
            
            a(1:2) = G_z.Denominator{1}(2:end);
            b(1:2) = G_z.Numerator{1}(2:end);
            K = dcgain(G_z);
        end

        function diff_eq(obj, a, b, u, y)
            y_mod = zeros(size(y));
            y_mod(1:obj.delay+2) = y(1:obj.delay+2);
        
            if strcmp(obj.mode, 'ARX')
                for k = obj.delay+3:obj.kk
                    y_mod(k) = - a*[y(k-1:-1:k-2)]' + b*[u(1, k-(obj.delay+1):-1:k-(obj.delay+2))]' + b*[u(2, k-1:-1:k-2)]';
                end
            else
                for k = obj.delay+3:obj.kk
                    y_mod(k) = - a*[y_mod(k-1:-1:k-2)]' + b*[u(1, k-(obj.delay+1):-1:k-(obj.delay+2))]' + b*[u(2, k-1:-1:k-2)]';
                end
            end
            
            E = sum((y - y_mod).^2)/(obj.kk);
            fprintf('Model %s \n', obj.mode);
            fprintf('E = %.3f \n', E);
        
            figure;
            hold on;
            stairs(0:obj.kk-1, y_mod, 'b-', 'LineWidth', 1.2);
            stairs(0:obj.kk-1, y, 'r-', 'LineWidth', 0.8);
            hold off;
            xlabel('k');
            ylabel('y(k)');
            title(sprintf('Model %s \n E = %.3f', obj.mode, E));
            legend('y_{mod}', 'y');
            grid on;
            % file_name = sprintf('../raport/pictures/arx_ucz.pdf');
            % exportgraphics (gcf, file_name);
        end

        function [y, y_L] = rk4(obj, u, kk)

            %% Alokacja pamięci
            y = zeros(1, kk);
            y_L = zeros(1, kk);
            
            %% Warunki początkowe
            V_1 = obj.V_10;
            V_2 = obj.V_20;
            V_1L = obj.V_10;
            V_2L = obj.V_20;
            
            %% Funkcje
            fun_1 = @(F_1, F_D, V_1) F_1 + F_D - obj.alpha_1 * (V_1/obj.A)^(1/2);
            fun_1L = @(F_1, F_D, V_1) F_1 + F_D - obj.alpha_1 * (obj.V_10/obj.A)^(1/2) - obj.alpha_1 / (2*obj.V_10^(1/2)*obj.A^(1/2)) * (V_1-obj.V_10);
            fun_2 = @(V_1, V_2) obj.alpha_1 * (V_1/obj.A)^(1/2) - obj.alpha_2 * (V_2/obj.C)^(1/4); 
            fun_2L = @(V_1, V_2) obj.alpha_1 * (obj.V_10/obj.A)^(1/2) - obj.alpha_2 * (obj.V_20/obj.C)^(1/4) + ...
                     obj.alpha_1 / (2*obj.V_10^(1/2)*obj.A^(1/2)) * (V_1-obj.V_10) - obj.alpha_2 / (4*obj.V_20^(3/4)*obj.C^(1/4)) * (V_2-obj.V_20);
            
            %% Wymuszenia
            F_1in = u(1,:) + obj.F_10;
            F_D = u(2,:) + obj.F_D0;
            
            %% RK4 - Runge-Kutta 4. rzędu
            for i = 2:kk
                if i <= obj.delay + 1
                    % RK4 dla nieliniowego układu
                    % Krok dla V_1
                    k1_V1 = obj.Tp * fun_1(obj.F_10, F_D(i), max(V_1, 0));
                    k2_V1 = obj.Tp * fun_1(obj.F_10, F_D(i), max(V_1 + 0.5 * k1_V1, 0));
                    k3_V1 = obj.Tp * fun_1(obj.F_10, F_D(i), max(V_1 + 0.5 * k2_V1, 0));
                    k4_V1 = obj.Tp * fun_1(obj.F_10, F_D(i), max( V_1 + k3_V1, 0));
                    V_1 = V_1 + (k1_V1 + 2*k2_V1 + 2*k3_V1 + k4_V1) / 6;
            
                    % Krok dla V_2
                    k1_V2 = obj.Tp * fun_2(max(V_1, 0), max(V_2, 0));
                    k2_V2 = obj.Tp * fun_2(max(V_1, 0), max(V_2 + 0.5 * k1_V2, 0));
                    k3_V2 = obj.Tp * fun_2(max(V_1, 0), max(V_2 + 0.5 * k2_V2, 0));
                    k4_V2 = obj.Tp * fun_2(max(V_1, 0), max(V_2 + k3_V2, 0));
                    V_2 = V_2 + (k1_V2 + 2*k2_V2 + 2*k3_V2 + k4_V2) / 6;
            
                    % RK4 dla układu liniowego
                    k1_V1L = obj.Tp * fun_1L(obj.F_10, F_D(i), V_1L);
                    k2_V1L = obj.Tp * fun_1L(obj.F_10, F_D(i), V_1L + 0.5 * k1_V1L);
                    k3_V1L = obj.Tp * fun_1L(obj.F_10, F_D(i), V_1L + 0.5 * k2_V1L);
                    k4_V1L = obj.Tp * fun_1L(obj.F_10, F_D(i), V_1L + k3_V1L);
                    V_1L = V_1L + (k1_V1L + 2*k2_V1L + 2*k3_V1L + k4_V1L) / 6;
            
                    k1_V2L = obj.Tp * fun_2L(V_1L, V_2L);
                    k2_V2L = obj.Tp * fun_2L(V_1L, V_2L + 0.5 * k1_V2L);
                    k3_V2L = obj.Tp * fun_2L(V_1L, V_2L + 0.5 * k2_V2L);
                    k4_V2L = obj.Tp * fun_2L(V_1L, V_2L + k3_V2L);
                    V_2L = V_2L + (k1_V2L + 2*k2_V2L + 2*k3_V2L + k4_V2L) / 6;
                else
                    % Uwzględnienie opóźnienia dla F_1in
                    delayed_idx = i - obj.delay;
            
                    % RK4 dla nieliniowego układu
                    k1_V1 = obj.Tp * fun_1(F_1in(delayed_idx), F_D(i), max(V_1, 0));
                    k2_V1 = obj.Tp * fun_1(F_1in(delayed_idx), F_D(i), max(V_1 + 0.5 * k1_V1, 0));
                    k3_V1 = obj.Tp * fun_1(F_1in(delayed_idx), F_D(i), max(V_1 + 0.5 * k2_V1, 0));
                    k4_V1 = obj.Tp * fun_1(F_1in(delayed_idx), F_D(i), max(V_1 + k3_V1, 0));
                    V_1 = V_1 + (k1_V1 + 2*k2_V1 + 2*k3_V1 + k4_V1) / 6;
            
                    k1_V2 = obj.Tp * fun_2(max(V_1, 0), max(V_2, 0));
                    k2_V2 = obj.Tp * fun_2(max(V_1, 0), max(V_2 + 0.5 * k1_V2, 0));
                    k3_V2 = obj.Tp * fun_2(max(V_1, 0), max(V_2 + 0.5 * k2_V2, 0));
                    k4_V2 = obj.Tp * fun_2(max(V_1, 0), max(V_2 + k3_V2, 0));
                    V_2 = V_2 + (k1_V2 + 2*k2_V2 + 2*k3_V2 + k4_V2) / 6;
            
                    % RK4 dla układu liniowego
                    k1_V1L = obj.Tp * fun_1L(F_1in(delayed_idx), F_D(i), V_1L);
                    k2_V1L = obj.Tp * fun_1L(F_1in(delayed_idx), F_D(i), V_1L + 0.5 * k1_V1L);
                    k3_V1L = obj.Tp * fun_1L(F_1in(delayed_idx), F_D(i), V_1L + 0.5 * k2_V1L);
                    k4_V1L = obj.Tp * fun_1L(F_1in(delayed_idx), F_D(i), V_1L + k3_V1L);
                    V_1L = V_1L + (k1_V1L + 2*k2_V1L + 2*k3_V1L + k4_V1L) / 6;
            
                    k1_V2L = obj.Tp * fun_2L(V_1L, V_2L);
                    k2_V2L = obj.Tp * fun_2L(V_1L, V_2L + 0.5 * k1_V2L);
                    k3_V2L = obj.Tp * fun_2L(V_1L, V_2L + 0.5 * k2_V2L);
                    k4_V2L = obj.Tp * fun_2L(V_1L, V_2L + k3_V2L);
                    V_2L = V_2L + (k1_V2L + 2*k2_V2L + 2*k3_V2L + k4_V2L) / 6;
                end
            
                V_2 = max(V_2, 0);
                V_2L = max(V_2L, 0);

                h_2e = sqrt(V_2 / obj.C);
                h_2eL = sqrt(obj.V_20 / obj.C) + 1/(2*sqrt(obj.V_20 * obj.C)) * (V_2L - obj.V_20);
            
                % Sprowadzenie wartości do punktu pracy
                y(i) = h_2e - obj.h_20;
                y_L(i) = h_2eL - obj.h_20;
            end
        end


        function show_sopdt(obj, y, G_z)
            %% Prezentacja wyników SOPDT
            figure;
            plot(0:obj.Tp:(obj.kk-1)*obj.Tp, y, 'r-','LineWidth',2);
            hold on;
            step(G_z);
            grid on;
            legend('h', 'h_{mod}', 'Location', 'northwest');
            % file_name = sprintf('../raport/pictures/model_sopdt.pdf');
            % exportgraphics (gcf, file_name);
        end

        function show_rk4(obj, u, y, y_L, kk)
            figure;
            stairs(0:kk-1, u(1,:), 'r-', 'LineWidth', 1.2);
            hold on;
            stairs(0:kk-1, u(2,:), 'b--', 'LineWidth', 1.2);
            hold off;
            xlabel('k');
            ylabel('u(k)');
            title('Wartości sygnałów sterujących u(k)');
            legend('F_1(k)', 'F_D(k)');
            grid on;

            %% Prezentacja wyników
            figure;
            hold on;
            plot(0:obj.Tp:(kk-1)*obj.Tp, y_L, 'b.','LineWidth',2);
            plot(0:obj.Tp:(kk-1)*obj.Tp, y, 'r-','LineWidth',2);
            hold off;
            title('Wysokość słupa cieczy w zbiorniku 2. - h_2(t)');
            legend('h_{2lin}', 'h_2');
            xlabel('t [s]');
            ylabel('h [cm]');
            grid on;
        end
    end
end