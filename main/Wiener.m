classdef Wiener < handle
    properties
        Y_min = -27;
        Y_max = 27;

        linear_fis
        nonlinear_fis
    end

    methods

        function obj = Wiener()
        end

        function linearFuzzy(obj)
            Y_center = linspace(obj.Y_min, obj.Y_max, 5);
            obj.linear_fis = sugfis('Name', 'Linear_Wiener', 'Type', 'sugeno');
            obj.linear_fis = addInput(obj.linear_fis, [obj.Y_min obj.Y_max], 'Name', 'Y_linear');
            
            % Definiowanie funkcji przynależności (gaussmf)
            for i = 1:length(Y_center) 
                obj.linear_fis = addMF(obj.linear_fis, 'Y_linear', 'gaussmf', [10, Y_center(i)]);
            end
            
            % Definiowanie wyjścia i początkowych następników (a_i * y + b_i)
            obj.linear_fis = addOutput(obj.linear_fis, [obj.Y_min obj.Y_max], 'Name', 'Y_fuzzy');
            
            % Początkowe współczynniki (a_i, b_i)
            a_param = [0.77 0.92 1 1.05 1.24]; % Można dostroić
            b_param = [0 0 0 0 0];
            
            % Dodanie reguł TS w postaci liniowej
            for i = 1:length(Y_center)
                obj.linear_fis = addMF(obj.linear_fis, 'Y_fuzzy', 'linear', [a_param(i), b_param(i)]);
            end
            
            % Reguły Takagi-Sugeno: [inputMF, outputMF, weight]
            ruleList = [1 1 1 1;
                        2 2 1 1;
                        3 3 1 1;
                        4 4 1 1;
                        5 5 1 1];
            
            % Dodanie reguł do systemu
            obj.linear_fis = addRule(obj.linear_fis, ruleList);
        end

        function nonlinearFuzzy(obj)
            Y_center = linspace(obj.Y_min, obj.Y_max, 3);
            obj.nonlinear_fis = sugfis('Name', 'Nonlinear_Wiener', 'Type', 'sugeno');
            obj.nonlinear_fis = addInput(obj.nonlinear_fis, [obj.Y_min obj.Y_max], 'Name', 'Y_linear');
            
            % Definiowanie funkcji przynależności (gaussmf)
            for i = 1:length(Y_center) 
                obj.nonlinear_fis = addMF(obj.nonlinear_fis, 'Y_linear', 'gaussmf', [15, Y_center(i)]);
            end
            
            % Definiowanie wyjścia i początkowych następników (a_i * y + b_i)
            obj.nonlinear_fis = addOutput(obj.nonlinear_fis, [obj.Y_min obj.Y_max], 'Name', 'Y_fuzzy');

            % Początkowe współczynniki (a_i, b_i)
            a_param = [44.1 47.2 50.7]; % Można dostroić
            b_param = [60 45 45];
                        
            % Dodanie reguł TS w postaci liniowej
            for i = 1:length(Y_center)
                obj.nonlinear_fis = addMF(obj.nonlinear_fis, 'Y_fuzzy', 'linear', [a_param(i), b_param(i)]);
            end
            
            % Reguły Takagi-Sugeno: [inputMF, outputMF, weight]
            ruleList = [1 1 1 1;
                        2 2 1 1;
                        3 3 1 1];
            
            % Dodanie reguł do systemu
            obj.nonlinear_fis = addRule(obj.nonlinear_fis, ruleList);
        end

        function testLinearModel(obj, U, a, b, delay, kk, Tp, rk4)
            [Y_real, Y_lin] = rk4(U, kk); % Symulacja rzeczywistego układu

            Y_out = zeros(1, kk);
            y_out = zeros(1, kk);
            for k = delay+3:kk
                y_out(k) = - a*[y_out(k-1:-1:k-2)]' + b*[U(1, k-(delay+1):-1:k-(delay+2))]' + b*[U(2, k-1:-1:k-2)]';
                Y_out(k) = evalfis(obj.linear_fis, y_out(k));
            end

            t = 0:Tp:(kk-1)*Tp;
            figure;
            plot(t, Y_real, 'b', t, Y_lin, 'g', t, Y_out, 'r');
            legend('RK4', 'RK4 liniowy', 'Wiener (liniowy TS)');
            title('Porównanie wyjścia układu rzeczywistego i modelu');
            ylabel('y [cm]');
            xlabel('t [s]');
            grid on;

            E_lin = sum((Y_real - Y_lin).^2);
            E_out = sum((Y_real - Y_out).^2);
            fprintf("\nWIENER LINEAR MODEL\n");
            fprintf("E_lin = %.3f\n", E_lin);
            fprintf("E_out = %.3f\n", E_out);
        end

        function testNonlinearModel(obj, U, a, b, delay, kk, Tp, rk4)
            [Y_real, Y_lin] = rk4(U, kk); % Symulacja rzeczywistego układu
            
            % Symulacja modelu Wienera
            y_out = zeros(1, kk);
            Y_out = zeros(1, kk);
            for k = delay+3:kk
                y_out(k) = - a*[y_out(k-1:-1:k-2)]' + b*[U(1, k-(delay+1):-1:k-(delay+2))]' + b*[U(2, k-1:-1:k-2)]';
                [~, degrees] = evalfis(obj.nonlinear_fis, y_out(k)); % Przepuszczenie przez model TS
                output = 0;
                for i = 1:length(obj.nonlinear_fis.Rules)
                    output = output + degrees(i) * (obj.nonlinear_fis.Outputs.MembershipFunctions(i).Parameters(1)*sinh(y_out(k)/obj.nonlinear_fis.Outputs.MembershipFunctions(i).Parameters(2)));
                end
                Y_out(k) = output / sum(degrees);
            end
            
            t = 0:Tp:(kk-1)*Tp;
            figure;
            plot(t, Y_real, 'b', t, Y_lin, 'g', t, Y_out, 'r');
            legend('RK4', 'RK4 liniowy', 'Wiener (optymalny TS)', 'Location', 'northwest');
            title('Porównanie wyjścia układu rzeczywistego i modelu');
            ylabel('y [cm]');
            xlabel('t [s]');
            grid on;
            
            E_lin = sum((Y_real - Y_lin).^2);
            E_out = sum((Y_real - Y_out).^2);
            fprintf("\nWIENER NONLINEAR MODEL\n");
            fprintf("E_lin = %.3f\n", E_lin);
            fprintf("E_out = %.3f\n", E_out);
        end


        function show_fuzzy(~, fis)
            figure;
            plotmf(fis, 'input', 1);
            title('Funkcje przynależności dla Y');
            xlabel('y');
            ylabel('$\mu(y)$', 'Interpreter', 'latex');
            grid on;
        end

    end
end
