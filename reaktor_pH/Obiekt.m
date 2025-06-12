classdef Obiekt < handle
    properties
        % Stałe
        A = 207;
        C_V = 8.75;
        pK_1 = 6.35;
        pK_2 = 10.25;
        W_a1 = 3*10^(-3);
        W_a2 = -3*10^(-2);
        W_a3 = -3.05*10^(-3);
        W_b1 = 0;
        W_b2 = 3*10^(-2);
        W_b3 = 5*10^(-5);

        % Opóźnienie
        tau = 30;
        % Okres próbkowania
        Tp = 10;
        % Próbki dyskretne
        kk = 2100;  % 1200

        % Punkt pracy
        % Q_1, Q_3 - wartości sterujące
        % Q_2 - zakłócenie
        % Q_4 - do obliczenia (grawitacyjny wypływ ze zbiornika)
        Q_10
        Q_20
        Q_30
        Q_40
        W_a40
        W_b40
        h_0
        pH_0

        % Współczynniki modelu
        delay
        dynamic_horizont = 120;
    end

    methods
        function obj = Obiekt()
            obj.delay = obj.tau/obj.Tp;
        end

        function linearization(obj, Q_10, Q_20, Q_30)

            obj.h_0 = ((Q_10+Q_20+Q_30) / (obj.C_V))^2;

            obj.Q_10 = Q_10;
            obj.Q_20 = Q_20;
            obj.Q_30 = Q_30;
            obj.Q_40 = obj.C_V * obj.h_0;

            obj.W_a40 = (obj.W_a1*Q_10 + obj.W_a2*Q_20 + obj.W_a3*Q_30) / (Q_10 + Q_20 + Q_30 );
            obj.W_b40 = (obj.W_b1*Q_10 + obj.W_b2*Q_20 + obj.W_b3*Q_30) / (Q_10 + Q_20 + Q_30 );

            obj.pH_0 = obj.pH_calc(obj.W_a40, obj.W_b40);
        end

        function [a, b, S] = fopdtModel(obj)

            u = [ones(1,obj.dynamic_horizont)
                zeros(1, obj.dynamic_horizont)
                zeros(1,obj.dynamic_horizont)];
            [s, ~] = obj.modifiedEuler(u, obj.dynamic_horizont);
            y = s(1,:) / s(1,end);

            t = (0:length(y)-1) * obj.Tp;
            
            % Funkcja celu dopasowania SOPDT (K=1, optymalizujemy tylko T1 i T2)
            cost_func = @(params) ...
                sum((step(tf(1, [params 1]), t) - y').^2);
            
            % Początkowe zgadywanie
            initial_guess = 200;  % [T1, T2]
            
            % Optymalizacja z fminsearch
            options = optimset('Display', 'final', 'TolX', 1e-6, 'TolFun', 1e-6);
            optimal_params = fminsearch(cost_func, initial_guess, options);
            
            % Transmitancja dopasowanego modelu
            G = tf(1, [optimal_params 1]);
            S.Q1 = step(G, t);
            G_z = c2d(G, obj.Tp, 'zoh');
            G_z.Variable = 'z^-1';
            a = G_z.Denominator{1}(2);
            b = G_z.Numerator{1}(2);

            S.Q2 = S.Q1;
            S.Q3 = S.Q1;

            % Rysowanie wykresu
            figure;
            plot(t, y, 'b', 'LineWidth', 2); hold on;
            plot(t, S.Q1, 'r--', 'LineWidth', 2);
            legend('Odpowiedź rzeczywista', 'Model SOPDT');
            xlabel('Czas');
            ylabel('Odpowiedź skokowa');
            title(sprintf('Dopasowanie SOPDT: T1 = %.2f', optimal_params));
            grid on;
        end

        function [a, b, S] = tfestModel(obj)

            u = [ones(1,obj.dynamic_horizont)
                zeros(1,obj.dynamic_horizont)
                zeros(1,obj.dynamic_horizont)];
            [s, ~] = obj.modifiedEuler(u, obj.dynamic_horizont);
            y.Q1 = s(2,:) / s(2,end);
            model.Q1 = tfest(iddata(y.Q1', u(1,:)', obj.Tp), 2, 1);

            u = [zeros(1,obj.dynamic_horizont)
                ones(1,obj.dynamic_horizont)
                zeros(1,obj.dynamic_horizont)];
            [s, ~] = obj.modifiedEuler(u, obj.dynamic_horizont);
            y.Q2 = s(2,:) / s(2,end);
            model.Q2 = tfest(iddata(y.Q2', u(2,:)', obj.Tp), 2, 1);

            u = [zeros(1,obj.dynamic_horizont)
                zeros(1,obj.dynamic_horizont)
                ones(1,obj.dynamic_horizont)];
            [s, ~] = obj.modifiedEuler(u, obj.dynamic_horizont);
            y.Q3 = s(2,:) / s(2,end);
            model.Q3 = tfest(iddata(y.Q3', u(3,:)', obj.Tp), 2, 1);

            t = (0:length(y.Q1)-1) * obj.Tp;

            G = [tf(model.Q1), tf(model.Q2), tf(model.Q3)];
            G.InputName = {'u1', 'u2', 'u3'};
            G.OutputName = {'y1'};
            [Y, ~] = step(G, t);
            S.Q1 = Y(:,1,1);  % odpowiedź na skok w wejściu 1
            S.Q2 = Y(:,1,2);  % odpowiedź na skok w wejściu 2
            S.Q3 = Y(:,1,3);  % odpowiedź na skok w wejściu 3

            G_z = c2d(G, obj.Tp, 'zoh');
            G_z.Variable = 'z^-1';
            a.Q1 = G_z.Denominator{1}(2:3);
            b.Q1 = G_z.Numerator{1}(2:3);
            a.Q2 = G_z.Denominator{2}(2:3);
            b.Q2 = G_z.Numerator{2}(2:3);
            a.Q3 = G_z.Denominator{3}(2:3);
            b.Q3 = G_z.Numerator{3}(2:3);

            % Rysowanie wykresu
            figure;
            subplot(1,3,1);
            plot(t, y.Q1, 'b', 'LineWidth', 2); 
            hold on;
            plot(t, S.Q1, 'r--', 'LineWidth', 2);
            legend('Odpowiedź rzeczywista', 'Model', 'Location', 'best');
            xlabel('Czas');
            ylabel('Odpowiedź skokowa');
            title(sprintf('Odpowiedź układu na wymuszenie jednostkowe Q1'));
            grid on;

            subplot(1,3,2);
            plot(t, y.Q2, 'b', 'LineWidth', 2); 
            hold on;
            plot(t, S.Q2, 'r--', 'LineWidth', 2);
            legend('Odpowiedź rzeczywista', 'Model', 'Location', 'best');
            xlabel('Czas');
            ylabel('Odpowiedź skokowa');
            title(sprintf('Odpowiedź układu na wymuszenie jednostkowe Q2'));
            grid on;

            subplot(1,3,3);
            plot(t, y.Q3, 'b', 'LineWidth', 2); 
            hold on;
            plot(t, S.Q3, 'r--', 'LineWidth', 2);
            legend('Odpowiedź rzeczywista', 'Model', 'Location', 'best');
            xlabel('Czas');
            ylabel('Odpowiedź skokowa');
            title(sprintf('Odpowiedź układu na wymuszenie jednostkowe Q3'));
            grid on;
        end

        function [y, y_L, E_h, E_pH] = modifiedEuler(obj, u, kk)
            %% Alokacja pamięci
            y = zeros(2, kk);
            y_L = zeros(2, kk);
            W_a4 = zeros(1, kk);
            W_b4 = zeros(1, kk);
            W_a4e = zeros(1, kk);
            W_b4e = zeros(1, kk);
            W_a4L = zeros(1, kk);
            W_b4L = zeros(1, kk);
            W_a4eL = zeros(1, kk);
            W_b4eL = zeros(1, kk);

            %% Warunki początkowe
            W_a4(1) = obj.W_a40;
            W_b4(1) = obj.W_b40;
            W_a4e(1) = obj.W_a40;
            W_b4e(1) = obj.W_b40;

            h = obj.h_0;
            h_e = obj.h_0;

            W_a4L(1) = obj.W_a40;
            W_b4L(1) = obj.W_b40;
            W_a4eL(1) = obj.W_a40;
            W_b4eL(1) = obj.W_b40;

            h_L = obj.h_0;
            h_eL = obj.h_0;

            %% Funkcje
            fun_1 = @(Q_1, Q_2, Q_3, h) (Q_1 + Q_2 + Q_3 - obj.C_V*sqrt(h)) / obj.A;
            fun_1L = @(Q_1, Q_2, Q_3, h) (Q_1 + Q_2 + Q_3 - obj.C_V*sqrt(obj.h_0)) / obj.A - obj.C_V/(2*sqrt(obj.h_0)*obj.A) * (h - obj.h_0);
            fun_2 = @(Q_1, Q_2, Q_3, h, W_a4) ((obj.W_a1-W_a4)*Q_1 + (obj.W_a2-W_a4)*Q_2 + (obj.W_a3-W_a4)*Q_3) / (obj.A * h);
            fun_2L = @(Q_1, Q_2, Q_3, h, W_a4) ((obj.W_a1-obj.W_a40)*Q_1 + (obj.W_a2-obj.W_a40)*Q_2 + (obj.W_a3-obj.W_a40)*Q_3) / (obj.A * obj.h_0) ...
                - (obj.Q_10 + obj.Q_20 + obj.Q_30) / (obj.A*obj.h_0) * (W_a4 - obj.W_a40) ...
                + ((obj.W_a1-obj.W_a40)*obj.Q_10 + (obj.W_a2-obj.W_a40)*obj.Q_20 + (obj.W_a3-obj.W_a40)*obj.Q_30) / (obj.A * obj.h_0^2) * (h - obj.h_0);
            fun_3 = @(Q_1, Q_2, Q_3, h, W_b4) ((obj.W_b1-W_b4)*Q_1 + (obj.W_b2-W_b4)*Q_2 + (obj.W_b3-W_b4)*Q_3) / (obj.A * h);
            fun_3L = @(Q_1, Q_2, Q_3, h, W_b4) ((obj.W_b1-obj.W_b40)*Q_1 + (obj.W_b2-obj.W_b40)*Q_2 + (obj.W_b3-obj.W_b40)*Q_3) / (obj.A * obj.h_0) ...
                - (obj.Q_10 + obj.Q_20 + obj.Q_30) / (obj.A*obj.h_0) * (W_b4 - obj.W_b40) ...
                + ((obj.W_b1-obj.W_b40)*obj.Q_10 + (obj.W_b2-obj.W_b40)*obj.Q_20 + (obj.W_b3-obj.W_b40)*obj.Q_30) / (obj.A * obj.h_0^2) * (h - obj.h_0);
            
            %% Wymuszenia
            Q_1 = u(1,:) + obj.Q_10;
            Q_2 = u(2,:) + obj.Q_20;
            Q_3 = u(3,:) + obj.Q_30;
            
            %% Modified Euler
            for k = 2:kk
                h = h + obj.Tp * fun_1(Q_1(k), Q_2(k), Q_3(k), h);
                W_a4(k) = W_a4(k-1) + obj.Tp * fun_2(Q_1(k), Q_2(k), Q_3(k), h, W_a4(k-1));
                W_b4(k) = W_b4(k-1) + obj.Tp * fun_3(Q_1(k), Q_2(k), Q_3(k), h, W_b4(k-1));

                h_e = h_e + 1/2 * obj.Tp * (fun_1(Q_1(k), Q_2(k), Q_3(k), h_e) + fun_1(Q_1(k), Q_2(k), Q_3(k), h));
                W_a4e(k) = W_a4e(k-1) + 1/2 * obj.Tp * (fun_2(Q_1(k), Q_2(k), Q_3(k), h_e, W_a4e(k-1)) + fun_2(Q_1(k), Q_2(k), Q_3(k), h, W_a4(k)));
                W_b4e(k) = W_b4e(k-1) + 1/2 * obj.Tp * (fun_3(Q_1(k), Q_2(k), Q_3(k), h_e, W_b4e(k-1)) + fun_3(Q_1(k), Q_2(k), Q_3(k), h, W_b4(k)));
                
                % if (k <= obj.delay)
                %     pH = obj.pH_calc(obj.W_a40, obj.W_b40);
                % else
                %     pH = obj.pH_calc(W_a4e(k - obj.delay), W_b4e(k - obj.delay));
                % end
                pH = obj.pH_calc(W_a4e(k), W_b4e(k));
                % Sprowadzenie wartości do punktu pracy
                y(1, k) = h_e - obj.h_0;
                y(2, k) = pH - obj.pH_0;
            end

            for k = 2:kk
                h_L = h_L + obj.Tp * fun_1L(Q_1(k), Q_2(k), Q_3(k), h_L);
                W_a4L(k) = W_a4L(k-1) + obj.Tp * fun_2L(Q_1(k), Q_2(k), Q_3(k), h_L, W_a4L(k-1));
                W_b4L(k) = W_b4L(k-1) + obj.Tp * fun_3L(Q_1(k), Q_2(k), Q_3(k), h_L, W_b4L(k-1));

                h_eL = h_eL + 1/2 * obj.Tp * (fun_1L(Q_1(k), Q_2(k), Q_3(k), h_eL) + fun_1L(Q_1(k), Q_2(k), Q_3(k), h_L));
                W_a4eL(k) = W_a4eL(k-1) + 1/2 * obj.Tp * (fun_2L(Q_1(k), Q_2(k), Q_3(k), h_eL, W_a4eL(k-1)) + fun_2L(Q_1(k), Q_2(k), Q_3(k), h_L, W_a4L(k)));
                W_b4eL(k) = W_b4eL(k-1) + 1/2 * obj.Tp * (fun_3L(Q_1(k), Q_2(k), Q_3(k), h_eL, W_b4eL(k-1)) + fun_3L(Q_1(k), Q_2(k), Q_3(k), h_L, W_b4L(k)));
                
                % if (k <= obj.delay)
                %     pH = obj.pH_calc(obj.W_a40, obj.W_b40);
                % else
                %     pH = obj.pH_calc(W_a4eL(k - obj.delay), W_b4eL(k - obj.delay));
                % end
                pH = obj.pH_calc(W_a4eL(k), W_b4eL(k));
                % Sprowadzenie wartości do punktu pracy
                y_L(1, k) = h_eL - obj.h_0;
                y_L(2, k) = pH - obj.pH_0;
            end

            E_h = sum((y(1,:)-y_L(1,:)).^2) / kk;
            E_pH = sum((y(2,:)-y_L(2,:)).^2) / kk;
        end

        function [pH] = pH_calc(obj, W_a4, W_b4)
            a_0 = -10^obj.pK_1;
            a_1 = W_a4 * 10^obj.pK_1-1;
            a_2 = 10^(obj.pK_1-14) + W_a4 + W_b4 - 10^-obj.pK_2;
            a_3 = 10^-14 + W_a4*10^-obj.pK_2 + 2*W_b4*10^-obj.pK_2;
            a_4  = 10^-(obj.pK_2+14);

            p_kand=roots([a_4 a_3 a_2 a_1 a_0]);
            p=0;
            for j=1:4
                if isreal(p_kand(j))&&(p<p_kand(j))
                    p=p_kand(j);
                end
            end
            pH=log10(p);
        end

        function show_staticCharacteristic(obj, Q2)
            U_min = -15;
            U_max = 15;

            U = [linspace(U_min, U_max, 100);
                linspace(U_min, U_max, 100)];
            [Q1_grid, Q3_grid] = meshgrid(U(1,:), U(2,:));  % kombinacje
            
            h = zeros(100,100);
            pH = zeros(100,100);
            
            % Petla po siatce sterowań
            for i = 1:length(U)
                for j = 1:length(U)
                    u = [ones(1, 200)*U(1,i);
                         ones(1, 200)*Q2;
                         ones(1, 200)*U(2,j)];
                
                    [y, ~, ~] = obj.modifiedEuler(u, 200);
                    h(i,j) =  y(1, end);
                    pH(i,j) = y(2, end); 
                end
            end
            
            % Rysuj 3D wykres
            figure;
            surf(Q1_grid, Q3_grid, pH);
            xlabel('Q_1 [ml/s]');
            ylabel('Q_3 [ml/s]');
            zlabel('pH');
            title('Wpływ dopływów Q_1 oraz Q_3 na stężenie substancji pH');
            shading interp;
            colorbar;
            view(-45, 30);
            
            figure;
            surf(Q1_grid, Q3_grid, h);
            xlabel('Q_1 [ml/s]');
            ylabel('Q_3 [ml/s]');
            zlabel('h [cm]');
            title('Wpływ dopływów Q_1 oraz Q_3 na wysokość słupa cieczy h');
            shading interp;
            colorbar;
        end
    end
end