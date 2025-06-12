classdef Obiekt < handle
    properties
        % Stałe
        A = 540;
        C = 0.85;
        alpha_1 = 26;
        alpha_2 = 20;
        
        % Opóźnienie
        tau = 100;
        % Okres próbkowania
        Tp = 20;
        % Próbki dyskretne
        kk = 2000;

        % Punkt pracy
        F_10;
        F_D0;
        
        % Do obliczenia
        h_10
        h_20
        V_10
        V_20

        % Współczynniki modelu
        delay
        dynamic_horizont = 150
    end

    methods
        function obj = Obiekt()
            obj.delay = obj.tau/obj.Tp;
        end

        function linearization(obj, F_10, F_D0)
            obj.F_10 = F_10;
            obj.F_D0 = F_D0;

            obj.h_10 = ((F_10+F_D0)/obj.alpha_1)^2;
            obj.h_20 = ((F_10+F_D0)/obj.alpha_2)^2;
            obj.V_10 = obj.A * obj.h_10;
            obj.V_20 = obj.C * obj.h_20^2;
        end

        function [a, b, S, S_D] = sopdt(obj)
            u = [ones(1,obj.dynamic_horizont)
                zeros(1,obj.dynamic_horizont)];
            [s, ~] = obj.modifiedEuler(u, obj.dynamic_horizont);
            y = s / s(end);

            t = (0:length(y)-1) * obj.Tp;
            
            % Funkcja celu dopasowania SOPDT (K=1, optymalizujemy tylko T1 i T2)
            cost_func = @(params) ...
                sum((step(tf(1, conv([params(1) 1], [params(2) 1]), 'InputDelay', obj.tau), t) - y').^2);
            
            % Początkowe zgadywanie
            initial_guess = [200, 20];  % [T1, T2]
            
            % Optymalizacja z fminsearch
            options = optimset('Display', 'final', 'TolX', 1e-6, 'TolFun', 1e-6);
            optimal_params = fminsearch(cost_func, initial_guess, options);
            
            % Transmitancja dopasowanego modelu
            G = tf(1, conv([optimal_params(1) 1], [optimal_params(2) 1]));
            S = step(G, t);
            G = tf(G.Numerator, G.Denominator, 'InputDelay', obj.tau);
            S_D = step(G, t);
            G_z = c2d(G, obj.Tp, 'zoh');
            G_z.Variable = 'z^-1';
            a = G_z.Denominator{1}(2:3);
            b = G_z.Numerator{1}(2:3);

            % % Rysowanie wykresu
            % figure;
            % plot(t, y, 'b', 'LineWidth', 2); hold on;
            % plot(t, S, 'r--', 'LineWidth', 2);
            % legend('Odpowiedź rzeczywista', 'Model SOPDT');
            % xlabel('Czas');
            % ylabel('Odpowiedź skokowa');
            % title(sprintf('Dopasowanie SOPDT: T1 = %.2f, T2 = %.2f', optimal_params(1), optimal_params(2)));
            % grid on;
        end

        function [y, y_L, E] = modifiedEuler(obj, u, kk)

            %% Alokacja pamięci
            y = zeros(1, kk);
            y_L = zeros(1, kk);

            %% Warunki początkowe
            V_1 = obj.V_10;
            V_2 = obj.V_20;
            V_1L = obj.V_10;
            V_2L = obj.V_20;
            V_1e = obj.V_10;
            V_2e = obj.V_20;
            V_1eL = obj.V_10;
            V_2eL = obj.V_20;
            y(1) = 0;
            y_L(1) = 0;

            %% Funkcje
            fun_1 = @(F_1, F_D, V_1) F_1 + F_D - obj.alpha_1 * (V_1/obj.A)^(1/2);
            fun_1L = @(F_1, F_D, V_1) F_1 + F_D - obj.alpha_1 * (obj.V_10/obj.A)^(1/2) - obj.alpha_1 / (2*obj.V_10^(1/2)*obj.A^(1/2)) * (V_1-obj.V_10);
            fun_2 = @(V_1, V_2) obj.alpha_1 * (V_1/obj.A)^(1/2) - obj.alpha_2 * (V_2/obj.C)^(1/4); 
            fun_2L = @(V_1, V_2) obj.alpha_1 * (obj.V_10/obj.A)^(1/2) - obj.alpha_2 * (obj.V_20/obj.C)^(1/4) + ...
                     obj.alpha_1 / (2*obj.V_10^(1/2)*obj.A^(1/2)) * (V_1-obj.V_10) - obj.alpha_2 / (4*obj.V_20^(3/4)*obj.C^(1/4)) * (V_2-obj.V_20);

            %% Wymuszenia
            F_1in = u(1,:) + obj.F_10;
            F_D = u(2,:) + obj.F_D0;

            %% Modified Euler
            for i = 2:kk
                if i <= obj.delay+1
                    V_1 = V_1 + obj.Tp * fun_1(obj.F_10, F_D(i), V_1);
                    V_2 = V_2 + obj.Tp * fun_2(V_1, V_2);

                    V_1e = V_1e + 1/2 * obj.Tp * (fun_1(obj.F_10, F_D(i), V_1e) + fun_1(obj.F_10, F_D(i), V_1));
                    V_2e = V_2e + 1/2 * obj.Tp * (fun_2(V_1e, V_2e) + fun_2(V_1, V_2));

                    V_1L = V_1L + obj.Tp * fun_1L(obj.F_10, F_D(i), V_1L);
                    V_2L = V_2L + obj.Tp * fun_2L(V_1L, V_2L);

                    V_1eL = V_1eL + 1/2 * obj.Tp * (fun_1L(obj.F_10, F_D(i), V_1eL) + fun_1L(obj.F_10, F_D(i), V_1L));
                    V_2eL = V_2eL + 1/2 * obj.Tp * (fun_2L(V_1eL, V_2eL) + fun_2L(V_1L, V_2L));
                else            
                    V_1 = V_1 + obj.Tp * fun_1(F_1in(i - (obj.delay+1)), F_D(i), V_1);
                    V_2 = V_2 + obj.Tp * fun_2(V_1, V_2);

                    V_1e = V_1e + 1/2 * obj.Tp * (fun_1(F_1in(i - (obj.delay+1)), F_D(i), V_1e) + fun_1(F_1in(i - (obj.delay+1)), F_D(i), V_1));
                    V_2e = V_2e + 1/2 * obj.Tp * (fun_2(V_1e, V_2e) + fun_2(V_1, V_2));

                    V_1L = V_1L + obj.Tp * fun_1L(F_1in(i - (obj.delay+1)), F_D(i), V_1L);
                    V_2L = V_2L + obj.Tp * fun_2L(V_1L, V_2L);

                    V_1eL = V_1eL + 1/2 * obj.Tp * (fun_1L(F_1in(i - (obj.delay+1)), F_D(i), V_1eL) + fun_1L(F_1in(i - (obj.delay+1)), F_D(i), V_1L));
                    V_2eL = V_2eL + 1/2 * obj.Tp * (fun_2L(V_1eL, V_2eL) + fun_2L(V_1L, V_2L));
                end

                h_2 = sqrt(V_2e / obj.C);
                h_2L = sqrt(obj.V_20 / obj.C) + 1/(2*sqrt(obj.V_20 * obj.C)) * (V_2eL - obj.V_20);

                % Sprowadzenie wartości do punktu pracy
                y(i) = h_2 - obj.h_20;
                y_L(i) = h_2L - obj.h_20;
            end
            E = sum((y-y_L).^2) / kk;
        end

        function static_charakteristic(obj)
            F_1 = linspace(45, 135, obj.kk);
            F_D = linspace(15, 45, obj.kk);

            [F1_grid, FD_grid] = meshgrid(F_1, F_D);
        
            y = ((F1_grid + FD_grid) / obj.alpha_2).^2;
        
            figure;
            surf(F1_grid, FD_grid, y);
            xlabel('F_1');
            ylabel('F_D');
            zlabel('y');
            title('Charakterystyka statyczna');
            shading interp; % opcjonalne: wygładzenie powierzchni
        end
    end
end