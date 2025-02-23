classdef Obiekt < handle
    properties
        % Stałe
        a = 540;
        c = 0.85;
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

        function [s, s_disturbance, a, b] = linearization(obj, F_10, F_D0)
            obj.F_10 = F_10;
            obj.F_D0 = F_D0;

            obj.h_10 = ((F_10+F_D0)/obj.alpha_1)^2;
            obj.h_20 = ((F_10+F_D0)/obj.alpha_2)^2;
            obj.V_10 = obj.a * obj.h_10;
            obj.V_20 = obj.c * obj.h_20^2;

            A = [-obj.alpha_1/(2*obj.V_10^(1/2)*obj.a^(1/2)), 0;
                 obj.alpha_1/(2*sqrt(obj.V_10*obj.a)), - obj.alpha_2/(4*obj.V_20^(3/4)*obj.c^(1/4))];
            B = [1 1; 0 0];
            C = [0 1/(2*sqrt(obj.V_20*obj.c))];
            D = [0 0];

            sys = ss(A,B,C,D);
            G_s = tf(sys);
            G_s(1) = tf(G_s.Numerator(1), G_s.Denominator(1), 'InputDelay', obj.tau);
            G_z = c2d(G_s, obj.Tp, 'zoh');

            s = step(G_z(1), (0:obj.dynamic_horizont-1) * obj.Tp);
            s_disturbance = step(G_z(2), (0:obj.dynamic_horizont-1) * obj.Tp);
            a = G_z(1).Denominator{1}(2:end);
            b = G_z(1).Numerator{1}(2:end);
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
            y(1) = 0;
            y_L(1) = 0;
            
            %% Funkcje
            fun_1 = @(F_1, F_D, V_1) F_1 + F_D - obj.alpha_1 * (V_1/obj.a)^(1/2);
            fun_1L = @(F_1, F_D, V_1) F_1 + F_D - obj.alpha_1 * (obj.V_10/obj.a)^(1/2) - obj.alpha_1 / (2*obj.V_10^(1/2)*obj.a^(1/2)) * (V_1-obj.V_10);
            fun_2 = @(V_1, V_2) obj.alpha_1 * (V_1/obj.a)^(1/2) - obj.alpha_2 * (V_2/obj.c)^(1/4); 
            fun_2L = @(V_1, V_2) obj.alpha_1 * (obj.V_10/obj.a)^(1/2) - obj.alpha_2 * (obj.V_20/obj.c)^(1/4) + ...
                     obj.alpha_1 / (2*obj.V_10^(1/2)*obj.a^(1/2)) * (V_1-obj.V_10) - obj.alpha_2 / (4*obj.V_20^(3/4)*obj.c^(1/4)) * (V_2-obj.V_20);
            
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
            
                % V_2 = max(V_2, 0);
                % V_2L = max(V_2L, 0);

                h_2 = sqrt(V_2 / obj.c);
                h_2L = sqrt(obj.V_20 / obj.c) + 1/(2*sqrt(obj.V_20 * obj.c)) * (V_2L - obj.V_20);
            
                % Sprowadzenie wartości do punktu pracy
                y(i) = h_2 - obj.h_20;
                y_L(i) = h_2L - obj.h_20;
            end
        end

        function [fis_trained] = fuzzyfication(obj)
            U = zeros(2, obj.kk);
            U(1,:) =  repelem([0, -22.5, -45, 22.5, 45], 400);
            U(2,:) = repelem([0 -5, 0, 5], 500);
            
            [Y, ~] = obj.rk4(U, obj.kk);
            Y = Y';
            U = U';
            
            X = [Y(1:end-1), [zeros(obj.delay,1); U(1:end-(obj.delay+1), 1)], U(1:end-1, 2)];
            Y = Y(2:end);
            
            options = genfisOptions('GridPartition'); 
            options.NumMembershipFunctions = [2, 2, 2]; % Liczba funkcji przynależności
            options.InputMembershipFunctionType = 'gaussmf'; % Typ funkcji przynależności
            fis = genfis(X, Y, options);
            
            options = anfisOptions('InitialFIS', fis, 'EpochNumber', 100, 'DisplayErrorValues', false);
            fis_trained = anfis([X Y], options);
            % % Symulacja modelu na danych testowych
            % Y_pred = evalfis(fis_trained, X);
            
            % figure;
            % plotmf(fis_trained, 'input', 2);
            % title('Funkcje przynależności dla F_1');
            % grid on;

            % figure;
            % plotmf(fis_trained, 'input', 3);
            % title('Funkcje przynależności dla F_D');
            % grid on;
            
            % Wizualizacja wyników
            % figure;
            % plot(1:length(Y), Y, 'b', 'LineWidth', 1.5); hold on;
            % plot(1:length(Y_pred), Y_pred, 'r--', 'LineWidth', 1.5);
            % legend('Rzeczywiste wyjście', 'Model TS');
            % xlabel('Próbki');
            % ylabel('Wartość wyjściowa');
            % title('Porównanie rzeczywistego wyjścia z modelem Takagi-Sugeno');
            % grid on;
            
            % % Liczba reguł w FIS
            % numRules = length(fis.Rules);
            % 
            % % Wyświetlenie szczegółów następników dla każdej reguły
            % for i = 1:numRules
            %     fprintf('Reguła %d:\n', i);
            %     fprintf('   Opis: %s\n', fis_trained.Rules(i).Description);
            % 
            %     % Współczynniki następników są w właściwości Outputs.MembershipFunctions
            %     coeffs = fis_trained.Outputs.MembershipFunctions(i).Parameters;
            %     fprintf('   Następnik: y = ');
            %     fprintf('%f*x + ', coeffs(1:end-1)); % Współczynniki wejściowe
            %     fprintf('%f\n', coeffs(end)); % Wyraz wolny
            % end
        end

        function show_rk4(obj, u, y, y_L)
            figure;
            stairs(0:obj.kk-1, u(1,:), 'r-', 'LineWidth', 1.2);
            hold on;
            stairs(0:obj.kk-1, u(2,:), 'b--', 'LineWidth', 1.2);
            hold off;
            xlabel('k');
            ylabel('u(k)');
            title('Wartości sygnałów sterujących u(k)');
            legend('F_1(k)', 'F_D(k)');
            grid on;

            %% Prezentacja wyników
            figure;
            hold on;
            plot(0:obj.Tp:(obj.kk-1)*obj.Tp, round(y_L,3), 'b.','LineWidth',2);
            plot(0:obj.Tp:(obj.kk-1)*obj.Tp, round(y,3), 'r-','LineWidth',2);
            hold off;
            title('Wysokość słupa cieczy w zbiorniku 2. - h_2(t)');
            legend('h_{2lin}', 'h_2');
            xlabel('t [s]');
            ylabel('h [cm]');
            grid on;
        end
    end
end