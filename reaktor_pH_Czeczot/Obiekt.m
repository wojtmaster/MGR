classdef Obiekt < handle
    properties
        % Stałe
        V = 3
        Ka = 1.8*10^(-5);
        Kw = 10^(-14);

        % Okres próbkowania
        Tp = 1;
        % Próbki dyskretne
        kk = 100;

        % Punkt pracy
        C_10;
        C_20;
        C_1in_0;
        C_2in_0;
        F_10;
        F_20;
        pH_10;
        pH_20;
        pH_0

        % Do obliczenia
        C_A0;
        C_B0;
        C_D0;
        T_0;
        T_C_0;

        % Współczynniki modelu
        dynamic_horizont = 150
    end

    methods
        function obj = Obiekt()

        end

        function [] = linearization(obj, C_10, C_20, C_1in_0, C_2in_0, F_10, F_20, pH_10, pH_20, pH_0)

            obj.C_10 = C_10;
            obj.C_20 = C_20;
            obj.C_1in_0 = C_1in_0;
            obj.C_2in_0 = C_2in_0;
            obj.F_10 = F_10;
            obj.F_20 = F_20;
            obj.pH_10 = pH_10;
            obj.pH_20 = pH_20;
            obj.pH_0 = pH_0;

            % C = obj.V/obj.F_0*(obj.d*obj.k_o2*obj.C_B0*exp(-obj.E_2/(obj.R*obj.T_0)));
            % fprintf("C_D0 = %.3f\n", C);

            % A = [-obj.alpha_1/(2*obj.V_10^(1/2)*obj.a^(1/2)), 0;
            %      obj.alpha_1/(2*sqrt(obj.V_10*obj.a)), - obj.alpha_2/(4*obj.V_20^(3/4)*obj.c^(1/4))];
            % B = [1 1; 0 0];
            % C = [0 1/(2*sqrt(obj.V_20*obj.c))];
            % D = [0 0];
            %
            % sys = ss(A,B,C,D);
            % G_s = tf(sys);
            % G_s(1) = tf(G_s.Numerator(1), G_s.Denominator(1), 'InputDelay', obj.tau);
            % G_z = c2d(G_s, obj.Tp, 'zoh');
            %
            % s = step(G_z(1), (0:obj.dynamic_horizont-1) * obj.Tp);
            % s_disturbance = step(G_z(2), (0:obj.dynamic_horizont-1) * obj.Tp);
            % a = G_z(1).Denominator{1}(2:end);
            % b = G_z(1).Numerator{1}(2:end);
        end

        function [y] = rk4(obj, u)

            %% Alokacja pamięci
            y = zeros(2, obj.kk);

            %% Warunki początkowe
            C_1 = obj.C_10;
            C_2 = obj.C_20;

            %% Wymuszenia
            C_1in = u(1,:) + obj.C_1in_0;
            F_1 = u(2,:) + obj.F_10;
            C_2in = u(3,:) + obj.C_2in_0;
            F_2 = u(4,:) + obj.F_20;

            %% Funkcje
            fun_1 = @(F_1, C_1in, C_1, F_2) F_1/obj.V * (C_1in - C_1) - F_2/obj.V *C_1;
            fun_2 = @(F_2, C_2in, C_2, F_1) F_2/obj.V * (C_2in - C_2) - F_1/obj.V *C_2;

            y(:, 1) = 0;

            %% RK4 - Runge-Kutta 4. rzędu
            for i = 2:obj.kk
                % RK4 dla nieliniowego układu
                k1_C1 = obj.Tp * fun_1(F_1(i), C_1in(i), max(C_1, 0), F_2(i));
                k2_C1 = obj.Tp * fun_1(F_1(i), C_1in(i), max(C_1 + 0.5 * k1_C1, 0), F_2(i));
                k3_C1 = obj.Tp * fun_1(F_1(i), C_1in(i), max(C_1 + 0.5 * k2_C1, 0), F_2(i));
                k4_C1 = obj.Tp * fun_1(F_1(i), C_1in(i), max(C_1 + k3_C1, 0), F_2(i));
                C_1 = C_1 + (k1_C1 + 2*k2_C1 + 2*k3_C1 + k4_C1) / 6;

                k1_C2 = obj.Tp * fun_2(F_2(i), C_2in(i), max(C_2, 0), F_1(i));
                k2_C2 = obj.Tp * fun_2(F_2(i), C_2in(i), max(C_2 + 0.5 * k1_C2, 0), F_1(i));
                k3_C2 = obj.Tp * fun_2(F_2(i), C_2in(i), max(C_2 + 0.5 * k2_C2, 0), F_1(i));
                k4_C2 = obj.Tp * fun_2(F_2(i), C_2in(i), max(C_2 + k3_C2, 0), F_1(i));
                C_2 = C_2 + (k1_C2 + 2*k2_C2 + 2*k3_C2 + k4_C2) / 6;

                % Sprowadzenie wartości do punktu pracy
                y(1, i) = C_1 - obj.C_10;
                y(2, i) = C_2 - obj.C_20;
            end
        end

        function [y] = modifiedEuler(obj, u)

            %% Alokacja pamięci
            syms x;
            y = zeros(3, obj.kk);

            %% Warunki początkowe
            C_1 = obj.C_10;
            C_2 = obj.C_20;
            C_1e = obj.C_10;
            C_2e = obj.C_20;

            %% Wymuszenia
            C_1in = u(1,:) + obj.C_1in_0;
            F_1 = u(2,:) + obj.F_10;
            C_2in = u(3,:) + obj.C_2in_0;
            F_2 = u(4,:) + obj.F_20;

            %% Funkcje
            fun_1 = @(F_1, C_1in, C_1, F_2) F_1/obj.V * (C_1in - C_1) - F_2/obj.V *C_1;
            fun_2 = @(F_2, C_2in, C_2, F_1) F_2/obj.V * (C_2in - C_2) - F_1/obj.V *C_2;

            y(:, 1) = 0;

            %% Modified Euler
            for i = 2:obj.kk
                C_1 = C_1 + obj.Tp * fun_1(F_1(i), C_1in(i), max(C_1, 0), F_2(i));
                C_2 = C_2 + obj.Tp * fun_2(F_2(i), C_2in(i), max(C_2, 0), F_1(i));

                C_1e = C_1e + 1/2 * obj.Tp * (fun_1(F_1(i), C_1in(i), max(C_1e, 0), F_2(i)) + fun_1(F_1(i), C_1in(i), max(C_1, 0), F_2(i)));
                C_2e = C_2e + 1/2 * obj.Tp * (fun_2(F_2(i), C_2in(i), max(C_2e, 0), F_1(i)) + fun_2(F_2(i), C_2in(i), max(C_2, 0), F_1(i)));

                eq = 10^(-3*x) + 10^(-2*x) * (obj.Ka + C_2e) + 10^(-x) * (obj.Ka * (C_2e - C_1e) - obj.Kw) - obj.Ka * obj.Kw == 0;

                % Rozwiązanie symboliczne
                sol = vpasolve(eq, x);

                % Zamiana na wartości numeryczne
                pH = double(sol);

                % Sprowadzenie wartości do punktu pracy
                y(1, i) = C_1e - obj.C_10;
                y(2, i) = C_2e - obj.C_20;
                y(3, i) = pH - obj.pH_0;
            end
        end

        function y = staticCharakteristic(obj, u, n)
            syms x;
            y = zeros(4, n);
            %% Funkcje
            fun_1 = @(F_1, C_1in, C_1, F_2) F_1/obj.V * (C_1in - C_1) - F_2/obj.V *C_1;
            fun_2 = @(F_2, C_2in, C_2, F_1) F_2/obj.V * (C_2in - C_2) - F_1/obj.V *C_2;

            for i = 1:4
                C_1 = obj.C_10;
                C_2 = obj.C_20;
                C_1e = obj.C_10;
                C_2e = obj.C_20;
                for j = 1:n
                    switch i
                        case 1
                            C_1in = u(1,j) + obj.C_1in_0;
                            F_1 = obj.F_10;
                            C_2in = obj.C_2in_0;
                            F_2 = obj.F_20;
                        case 2
                            C_1in = obj.C_1in_0;
                            F_1 =  u(2,j) + obj.F_10;
                            C_2in = obj.C_2in_0;
                            F_2 = obj.F_20;
                        case 3
                            C_1in = obj.C_1in_0;
                            F_1 =  obj.F_10;
                            C_2in = u(3,j) + obj.C_2in_0;
                            F_2 = obj.F_20;
                        case 4
                            C_1in = obj.C_1in_0;
                            F_1 =  obj.F_10;
                            C_2in = obj.C_2in_0;
                            F_2 = u(4,j) + obj.F_20;
                    end

                    C_1 = C_1 + obj.Tp * fun_1(F_1, C_1in, max(C_1, 0), F_2);
                    C_2 = C_2 + obj.Tp * fun_2(F_2, C_2in, max(C_2, 0), F_1);

                    C_1e = C_1e + 1/2 * obj.Tp * (fun_1(F_1, C_1in, max(C_1e, 0), F_2) + fun_1(F_1, C_1in, max(C_1, 0), F_2));
                    C_2e = C_2e + 1/2 * obj.Tp * (fun_2(F_2, C_2in, max(C_2e, 0), F_1) + fun_2(F_2, C_2in, max(C_2, 0), F_1));

                    eq = 10^(-3*x) + 10^(-2*x) * (obj.Ka + C_2e) + 10^(-x) * (obj.Ka * (C_2e - C_1e) - obj.Kw) - obj.Ka * obj.Kw == 0;
                    sol = vpasolve(eq, x);
                    pH = double(sol);

                    y(i, j) = pH;
                end
            end
        end

        % function [pH1, pH2, pH] = calculate_pH(obj, C1, C2, C1in, C2in, pH1, pH2, pH)
        %
        %     % Definicja układu równań
        %     equations = @(x) [
        %         10^(-3*x(1)) + 10^(-2*x(1)) * obj.Ka + 10^(-x(1)) * (-obj.Ka * C1in - obj.Kw) - obj.Ka * obj.Kw;
        %         10^(-3*x(2)) + 10^(-2*x(2)) * (obj.Ka + C2in) + 10^(-x(2)) * (obj.Ka * C2in - obj.Kw) - obj.Ka * obj.Kw;
        %         10^(-3*x(3)) + 10^(-2*x(3)) * (obj.Ka + C2) + 10^(-x(3)) * (obj.Ka * (C2 - C1) - obj.Kw) - obj.Ka * obj.Kw
        %     ];
        %
        %     % Wartości początkowe dla solvera
        %     x0 = [pH1, pH2, pH]; % Początkowe przybliżenia pH
        %     % x0 = [0, 0, 0]; % Początkowe przybliżenia pH
        %
        %     % Rozwiązywanie układu równań
        %     options = optimoptions('fsolve', ...
        %             'Display', 'off', ... % Wyświetla kolejne kroki iteracji
        %             'TolFun', 1e-12, ...   % Zmniejsza tolerancję funkcji (mniejsze błędy)
        %             'TolX', 1e-12, ...     % Zmniejsza tolerancję zmiennej (dokładniejsze `x`)
        %             'MaxIterations', 5000, ... % Zwiększa liczbę iteracji
        %             'Algorithm', 'trust-region-reflective'); % Alternatywny algorytm
        %     [sol, fval, exitflag]= fsolve(equations, x0, options);
        %     if (exitflag == -1)
        %         disp(['Exit flag: ', num2str(exitflag)]);
        %         disp(['Residuum: ', num2str(norm(fval))]);
        %     end
        %
        %     % Zapisanie wyników
        %     pH1 = sol(1);
        %     pH2 = sol(2);
        %     pH = sol(3);
        % end

        function [y] = calculate_pH(obj, pH, C1, C2)
            x(1) = pH;
            y(1) = 10^(-3*x(1)) + 10^(-2*x(1)) * (obj.Ka + C2) + 10^(-x(1)) * (obj.Ka * (C2 - C1) - obj.Kw) - obj.Ka * obj.Kw;
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

        function show_rk4(obj, u, y)
            figure;
            stairs(0:obj.kk-1, u(1,:), 'r-', 'LineWidth', 1.2);
            hold on;
            stairs(0:obj.kk-1, u(2,:), 'b--', 'LineWidth', 1.2);
            stairs(0:obj.kk-1, u(3,:), 'g--', 'LineWidth', 1.2);
            stairs(0:obj.kk-1, u(4,:), 'y--', 'LineWidth', 1.2);
            stairs(0:obj.kk-1, u(5,:), 'k--', 'LineWidth', 1.2);
            hold off;
            xlabel('k');
            ylabel('u(k)');
            title('Wartości sygnałów sterujących oraz zakłóceń');
            legend('C_Ain(k)', 'F(k)', 'F_C(k)', 'T_in(k)', 'T_Cin(k)');
            grid on;

            %% Prezentacja wyników
            % y_figure = figure;
            figure(y_figure);
            yyaxis left
            plot(0:obj.Tp:(obj.kk-1)*obj.Tp, y(1,:), 'b.','LineWidth',2);
            ylabel('C_D');

            yyaxis right
            plot(0:obj.Tp:(obj.kk-1)*obj.Tp, y(2,:), 'r-','LineWidth',2);
            ylabel('T');

            title('Wartości sygnałów wyjściowych y(k)');
            xlabel('t [s]');
            legend('C_D', 'T');
            grid on;
        end
    end
end