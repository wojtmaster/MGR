classdef Hammerstein < handle

    properties
        U_min = -45;
        U_max = 45;

        linear_fis
        nonlinear_fis
    end

    methods
        function obj = Hammerstein()
        end

        function linearFuzzy(obj)
            U_center = linspace(obj.U_min, obj.U_max, 5);
            obj.linear_fis = sugfis('Name', 'Linear_Hammerstein', 'Type', 'sugeno');
            obj.linear_fis = addInput(obj.linear_fis, [obj.U_min obj.U_max], 'Name', 'U_linear');
            
            % Definiowanie funkcji przynależności (gaussmf)
            for i = 1:length(U_center)
                obj.linear_fis = addMF(obj.linear_fis, 'U_linear', 'gaussmf', [12, U_center(i)]);
            end
            
            % Definiowanie wyjścia i początkowych następników (a_i * u + b_i)
            obj.linear_fis = addOutput(obj.linear_fis, [obj.U_min obj.U_max], 'Name', 'U_fuzzy');
            
            % Współczynniki (a_i, b_i) następników
            a_param = [0.7895 0.8982 1.0933 1.0489 1.2034];
            b_param = [0.0001 0.0002 0.0001 0 0.001];
            
            % Dodanie reguł TS w postaci liniowej
            for i = 1:length(U_center)
                obj.linear_fis = addMF(obj.linear_fis, 'U_fuzzy', 'linear', [a_param(i), b_param(i)]);
            end
            
            % Reguły Takagi-Sugeno: [inputMF, outputMF, weight]
            ruleList = [1 1 1 1;  % Reguła 1: wejście MF1 -> wyjście Out1
                        2 2 1 1;  % Reguła 2: wejście MF2 -> wyjście Out2
                        3 3 1 1;  % Reguła 3: wejście MF3 -> wyjście Out3
                        4 4 1 1;  % Reguła 4: wejście MF4 -> wyjście Out4
                        5 5 1 1]; % Reguła 5: wejście MF5 -> wyjście Out5
            
            % Dodanie reguł do systemu
            obj.linear_fis = addRule(obj.linear_fis, ruleList);
        end

        function nonlinearFuzzy(obj)
            U_center = linspace(obj.U_min, obj.U_max, 3);
            obj.nonlinear_fis = sugfis('Name', 'Nonlinear_Hammerstein', 'Type', 'sugeno');
            obj.nonlinear_fis = addInput(obj.nonlinear_fis, [obj.U_min obj.U_max], 'Name', 'U_linear');
            
            % Definiowanie funkcji przynależności (gaussmf)
            for i = 1:length(U_center)
                obj.nonlinear_fis = addMF(obj.nonlinear_fis, 'U_linear', 'gaussmf', [20, U_center(i)]);
            end
            
            % Definiowanie wyjścia i początkowych następników (a_i * u + b_i)
            obj.nonlinear_fis = addOutput(obj.nonlinear_fis, [obj.U_min obj.U_max], 'Name', 'U_fuzzy');

            % Początkowe współczynniki (a_i, b_i)
            % a_param = [60.5 50.7 83.5]; % Można dobrać inaczej
            % b_param = [80 50 75];
            a_param = [46.4941   36.2347  110.4841]; % Można dobrać inaczej
            b_param = [62.8862   36.2592   98.7690];
            % a_param = [77 69 65.2];
            % b_param = [100 67 60];
            
            % Dodanie reguł TS w postaci liniowej
            for i = 1:length(U_center)
                obj.nonlinear_fis = addMF(obj.nonlinear_fis, 'U_fuzzy', 'linear', [a_param(i), b_param(i)]);
            end
            
            % Reguły Takagi-Sugeno: [inputMF, outputMF, weight]
            ruleList = [1 1 1 1;  % Reguła 1: wejście MF1 -> wyjście Out1
                        2 2 1 1;  % Reguła 2: wejście MF2 -> wyjście Out2
                        3 3 1 1];  % Reguła 3: wejście MF3 -> wyjście Out3
            
            % Dodanie reguł do systemu
            obj.nonlinear_fis = addRule(obj.nonlinear_fis, ruleList);
        end

        function testLinearModel(obj, U, a, b, delay, kk, Tp, rk4, index)
            % U = [repelem((rand(1, kk/400) * 90 - 45), 400); repelem((rand(1, kk/250) * 10 - 5), 250)];
            [Y_real, Y_lin] = rk4(U, kk); % Symulacja rzeczywistego układu

            U_fuzzy = evalfis(obj.linear_fis, U(1,:)'); % Przepuszczenie przez model TS
            Y_out = zeros(1, kk);
            for k = delay+3:kk
                Y_out(k) = - a*[Y_out(k-1:-1:k-2)]' + b*[U_fuzzy(k-(delay+1):-1:k-(delay+2))] + b*[U(2, k-1:-1:k-2)]';
            end

            E_lin = sum((Y_real - Y_lin).^2) / kk;
            E_out = sum((Y_real - Y_out).^2) / kk;
            fprintf("\nHAMMERSTEIN LINEAR MODEL\n");
            fprintf("E_lin = %.3f\n", E_lin);
            fprintf("E_out = %.3f\n", E_out);

            t = 0:Tp:(kk-1)*Tp;
            figure;
            plot(t, Y_real, 'b', t, Y_lin, 'g', t, Y_out, 'r', 'LineWidth', 1.5);
            
            legend({'model nieliniowy', ...
                    ['model liniowy' blanks(17) sprintf('E = %.3f', E_lin)], ...
                    ['model Hammersteina' blanks(5) sprintf('E = %.3f', E_out)]}, ...
                    'Location', 'best');

            title('Porównanie wyjścia obiektu testowego i jego modeli');
            ylabel('y [cm]');
            xlabel('t [s]');
            grid on;

            % saveas(gcf, sprintf('D:/EiTI/MGR/raporty/raport_MGR/pictures/HammersteinLinearModel_%d.png', index));  % Zapisuje jako plik PNG
        end

        function testNonlinearModel(obj, U, a, b, delay, kk, Tp, rk4, index)
            U_fuzzy = zeros(1, kk);

            [Y_real, Y_lin] = rk4(U, kk); % Symulacja rzeczywistego układu
            Y_out = zeros(1, kk);

            for k = 1:kk
                [~, degrees] = evalfis(obj.nonlinear_fis, U(1, k));
                output = 0;
                for i = 1:length(obj.nonlinear_fis.Rules)
                    output = output + degrees(i) * (obj.nonlinear_fis.Outputs.MembershipFunctions(i).Parameters(1)*sinh(U(1,k)/obj.nonlinear_fis.Outputs.MembershipFunctions(i).Parameters(2)));
                end
                U_fuzzy(k) = output / sum(degrees);
                if k < delay+3
                    Y_out(k) = 0;
                else
                    Y_out(k) = - a*[Y_out(k-1:-1:k-2)]' + b*[U_fuzzy(k-(delay+1):-1:k-(delay+2))]' + b*[U(2, k-1:-1:k-2)]';
                end
            end

            E_lin = sum((Y_real - Y_lin).^2) / kk;
            E_out = sum((Y_real - Y_out).^2) / kk;
            fprintf("\nHAMMERSTEIN NONLINEAR MODEL\n");
            fprintf("E_lin = %.3f\n", E_lin);
            fprintf("E_out = %.3f\n", E_out);

            t = 0:Tp:(kk-1)*Tp;
            figure;
            plot(t, Y_real, 'b', t, Y_lin, 'g', t, Y_out, 'r', 'LineWidth', 1.5);
            legend({'model nieliniowy', ...
                    ['model liniowy' blanks(17) sprintf('E = %.3f', E_lin)], ...
                    ['model Hammersteina' blanks(5) sprintf('E = %.3f', E_out)]}, ...
                    'Location', 'best');
            title('Porównanie wyjścia obiektu testowego i jego modeli');
            ylabel('y [cm]');
            xlabel('t [s]');
            grid on;

            % saveas(gcf, sprintf('D:/EiTI/MGR/raporty/raport_MGR/pictures/HammersteinNonlinearModel_%d.png', index));  % Zapisuje jako plik PNG
        end

        function show_fuzzy(~, fis, s)
            figure;
            plotmf(fis, 'input', 1);
            title(sprintf('Funkcje przynależności dla u(k) - następniki %s', s));
            xlabel('u');
            ylabel('$\mu(u)$', 'Interpreter', 'latex');
            grid on;

            % saveas(gcf, sprintf('D:/EiTI/MGR/raporty/raport_MGR/pictures/HammersteinfuzzySets_%s.png', s));  % Zapisuje jako plik PNG
        end
    end
end