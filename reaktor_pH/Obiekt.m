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

            % a_1, a_2, a_3, b_1, b_2, b_3, G_z
            % A_matrix = [-obj.C_V/(2*sqrt(obj.h_0)*obj.A), 0, 0;
            %      ((obj.W_a1-obj.W_a40)*obj.Q_10 + (obj.W_a2-obj.W_a40)*obj.Q_20 + (obj.W_a3-obj.W_a40)*obj.Q_30) / (obj.A * obj.h_0^2), - (obj.Q_10 + obj.Q_20 + obj.Q_30) / (obj.A*obj.h_0), 0;
            %      ((obj.W_b1-obj.W_b40)*obj.Q_10 + (obj.W_b2-obj.W_b40)*obj.Q_20 + (obj.W_b3-obj.W_b40)*obj.Q_30) / (obj.A * obj.h_0^2), 0, - (obj.Q_10 + obj.Q_20 + obj.Q_30) / (obj.A*obj.h_0)];
            % B_matrix = [1 / obj.A, 1 / obj.A, 1 / obj.A; 
            %     (obj.W_a1-obj.W_a40) / (obj.A * obj.h_0), (obj.W_a2-obj.W_a40) / (obj.A * obj.h_0),  (obj.W_a3-obj.W_a40) / (obj.A * obj.h_0);
            %     (obj.W_b1-obj.W_b40) / (obj.A * obj.h_0),  (obj.W_b2-obj.W_b40) / (obj.A * obj.h_0), (obj.W_b3-obj.W_b40) / (obj.A * obj.h_0)];
            % C_matrix = eye(3);
            % D_matrix = zeros(3);
            % 
            % sys = ss(A_matrix, B_matrix, C_matrix, D_matrix);
            % G_s = tf(sys);
            % G_z = c2d(G_s, obj.Tp, 'zoh');
            % G_z = set(G_z, 'Variable', 'z^-1');  % poprawne ustawienie zmiennej
            % 
            % % s = step(G_z(1), (0:obj.dynamic_horizont-1) * obj.Tp);
            % a_1 = G_z(1,1).Denominator{1}(2:end);
            % a_2 = G_z(2,2).Denominator{1}(2:end);
            % a_3 = G_z(3,3).Denominator{1}(2:end);
            % b_1.Q_1 = G_z(1,1).Numerator{1}(2:end);
            % b_1.Q_2 = G_z(1,2).Numerator{1}(2:end);
            % b_1.Q_3 = G_z(1,3).Numerator{1}(2:end);
            % b_2.Q_1 = G_z(2,1).Numerator{1}(2:end);
            % b_2.Q_2 = G_z(2,2).Numerator{1}(2:end);
            % b_2.Q_3 = G_z(2,3).Numerator{1}(2:end);
            % b_3.Q_1 = G_z(3,1).Numerator{1}(2:end);
            % b_3.Q_2 = G_z(3,2).Numerator{1}(2:end);
            % b_3.Q_3 = G_z(3,3).Numerator{1}(2:end);
        end

        function [a, b, s] = mse(obj, data)
            if (strcmp(data, 'h'))
                u = [ones(1,obj.dynamic_horizont)
                    zeros(1,obj.dynamic_horizont)
                    zeros(1,obj.dynamic_horizont)];
                
                [s.Q1, ~] = obj.modifiedEuler(u, obj.dynamic_horizont);
                y.Q1 = s.Q1(1,:) / abs(s.Q1(1,end));
    
                Y = y.Q1(2:end)';
                M = [y.Q1(1:end-1)' u(1,1:end-1)'];
                w = M \ Y;
    
                a.Q1 = -w(1);
                b.Q1 = w(2);
    
                y_Q1 = zeros(size(y.Q1));
                y_Q1(1:2) = y.Q1(1:2);
                for k = 2:obj.dynamic_horizont
                    y_Q1(k) = -a.Q1*y_Q1(k-1) + b.Q1*u(1, k-1);
                end

                u = [zeros(1,obj.dynamic_horizont)
                    ones(1,obj.dynamic_horizont)
                    zeros(1,obj.dynamic_horizont)];
                [s.Q2, ~] = obj.modifiedEuler(u, obj.dynamic_horizont);
                y.Q2 = s.Q2(1,:);
    
                Y = y.Q2(2:end)';
                M = [y.Q2(1:end-1)' u(2,1:end-1)'];
                w = M \ Y;
    
                a.Q2 = -w(1);
                b.Q2 = w(2);
    
                y_Q2 = zeros(size(y.Q2));
                y_Q2(1:2) = y.Q2(1:2);
                for k = 2:obj.dynamic_horizont
                    y_Q2(k) = -a.Q2*y_Q2(k-1) + b.Q2*u(2, k-1);
                end

                u = [zeros(1,obj.dynamic_horizont)
                    zeros(1,obj.dynamic_horizont)
                    ones(1,obj.dynamic_horizont)];
                [s.Q3, ~] = obj.modifiedEuler(u, obj.dynamic_horizont);
                y.Q3 = s.Q3(1,:) / abs(s.Q3(1,end));
    
                Y = y.Q3(2:end)';
                M = [y.Q3(1:end-1)' u(3,1:end-1)'];
                w = M \ Y;
    
                a.Q3 = -w(1);
                b.Q3 = w(2);
    
                y_Q3 = zeros(size(y.Q3));
                y_Q3(1:2) = y.Q3(1:2);
                for k = 2:obj.dynamic_horizont
                    y_Q3(k) = -a.Q3*y_Q3(k-1) + b.Q3*u(3, k-1);
                end

            else
                u = [ones(1,obj.dynamic_horizont)
                    zeros(1,obj.dynamic_horizont)
                    zeros(1,obj.dynamic_horizont)];
                
                [s.Q1, ~] = obj.modifiedEuler(u, obj.dynamic_horizont);
                y.Q1 = s.Q1(2,:) / abs(s.Q1(2,end));
    
                Y = y.Q1(3:end)';
                M = [y.Q1(2:end-1)' y.Q1(1:end-2)' u(1,2:end-1)'];
                w = M \ Y;
    
                a.Q1 = -w(1:2)';
                b.Q1 = w(3)';
    
                y_Q1 = zeros(size(y.Q1));
                y_Q1(1:3) = y.Q1(1:3);
                for k = 4:obj.dynamic_horizont
                    y_Q1(k) = -a.Q1*y_Q1(k-1:-1:k-2)' + b.Q1*u(1, k-1);
                end

                u = [zeros(1,obj.dynamic_horizont)
                    ones(1,obj.dynamic_horizont)
                    zeros(1,obj.dynamic_horizont)];

                [s.Q2, ~] = obj.modifiedEuler(u, obj.dynamic_horizont);
                y.Q2 = s.Q2(2,:);

                Y = y.Q2(3:end)';
                M = [y.Q2(2:end-1)' y.Q2(1:end-2)' u(2,2:end-1)'];
                w = M \ Y;

                a.Q2 = -w(1:2)';
                b.Q2 = w(3);

                y_Q2 = zeros(size(y.Q2));
                y_Q2(1:3) = y.Q2(1:3);
                for k = 4:obj.dynamic_horizont
                    y_Q2(k) = -a.Q2*y_Q2(k-1:-1:k-2)' + b.Q2*u(2, k-1);
                end

                u = [zeros(1,obj.dynamic_horizont)
                    zeros(1,obj.dynamic_horizont)
                    ones(1,obj.dynamic_horizont)];

                [s.Q3, ~] = obj.modifiedEuler(u, obj.dynamic_horizont);
                y.Q3 = s.Q3(2,:) / abs(s.Q3(2,end));

                Y = y.Q3(4:end)';
                M = [y.Q3(3:end-1)' y.Q3(2:end-2)' y.Q3(1:end-3)' u(3,3:end-1)'];
                w = M \ Y;

                a.Q3 = -w(1:3)';
                b.Q3 = w(4)';

                y_Q3 = zeros(size(y.Q3));
                y_Q3(1:4) = y.Q3(1:4);
                for k = 5:obj.dynamic_horizont
                    y_Q3(k) = -a.Q3*y_Q3(k-1:-1:k-3)' + b.Q3*u(3, k-1);
                end

                % u = [repelem((rand(1, obj.kk/300) * 30 - 15), 300)
                %     zeros(1, obj.kk);
                %     repelem((rand(1, obj.kk/700) * 30 - 15), 700)];
                % 
                % [s.Q, ~] = obj.modifiedEuler(u, obj.kk);
                % y.Q = s.Q(2,:);
                % 
                % disp(size(y.Q));
                % disp(size(u));
                % 
                % Y = y.Q(4:end)';
                % M = [y.Q(3:end-1)' y.Q(2:end-2)' y.Q(1:end-3)' u(1, 3:end-1)' u(1, 2:end-2)' u(1, 1:end-3)' u(3, 3:end-1)' u(3, 2:end-2)' u(3, 1:end-3)'];
                % w = M \ Y;
                % 
                % a.Q = -w(1:3)';
                % b.Q = w(4:9)';
                % 
                % y_Q = zeros(size(y.Q));
                % y_Q(1:4) = y.Q(1:4);
                % for k = 5:obj.kk
                %     y_Q(k) = -[-1.9824 1.2739 -0.2690]*y_Q(k-1:-1:k-3)' + -0.02*u(1, k-1)' + 0.001*u(3, k-1)';
                % end
            end
            
            t = (0:length(u)-1) * obj.Tp;

            figure; 
            plot(t, y.Q1, 'b', t, y_Q1, 'g');
            title(sprintf('%s ( Q_1 )', data));
            legend('dane', 'symulacja');
            grid on;

            figure; 
            plot(t, y.Q2, 'm', t, y_Q2, 'c');
            title(sprintf('%s ( Q_2 )', data));
            legend('dane', 'symulacja');
            grid on;

            figure; 
            plot(t, y.Q3, 'y', t, y_Q3, 'k');
            title(sprintf('%s ( Q_3 )', data));
            legend('dane', 'symulacja');
            grid on;

            % figure; 
            % plot(t, y.Q, 'y', t, y_Q, 'k');
            % title(sprintf('%s ( Q_3 )', data));
            % legend('dane', 'symulacja');
            % grid on;
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