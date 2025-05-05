classdef Hammerstein < handle

    properties
        F1_min = -45;
        F1_max = 45;
        FD_min = -15;
        FD_max = 15;

        linear_fis
        nonlinear_fis
    end

    methods
        function obj = Hammerstein()
        end

        function linearFuzzy(obj)
            F1_center = linspace(-45, 45, 5);
            FD_center = linspace(-15, 15, 3);
            obj.linear_fis = sugfis('Name', 'Hammerstein', 'Type', 'sugeno');
            obj.linear_fis = addInput(obj.linear_fis, [obj.F1_min obj.F1_max], 'Name', 'F1');
            obj.linear_fis = addInput(obj.linear_fis, [obj.FD_min obj.FD_max], 'Name', 'FD');
            
            % Definiowanie funkcji przynależności (gaussmf)
            for i = 1:length(F1_center) 
                obj.linear_fis = addMF(obj.linear_fis, 'F1', 'gaussmf', [12, F1_center(i)]);
            end
            for i = 1:length(FD_center) 
                obj.linear_fis = addMF(obj.linear_fis, 'FD', 'gaussmf', [8, FD_center(i)]);
            end
            
            % Definiowanie wyjścia i początkowych następników (a_i * u + b_i)
            obj.linear_fis = addOutput(obj.linear_fis, [-27 45], 'Name', 'Y_fuzzy');
           
            a_param = [	0.3677	0.4958	0.7417	0.4476	0.6090	0.7706	0.5915	0.5585	0.7681	0.5531	0.5410	0.6824	0.5186	0.7147	0.6810];
            b_param = [	0.7801	0.8602	1.0268	0.5732	0.7117	0.7748	0.7151	0.5053	0.6096	0.6733	0.6520	0.5722	0.4357	0.5604	0.8445];
            c_param = [	1.4143	2.2959	1.9559	1.8199	0.4147	1.6517	0.6147	0.3596	-0.6468	1.8442	0.9703	1.6752	3.1648	2.5298	2.8351];

            % a_param = [0.5273    0.5902    0.5807    0.6200    0.6104    0.6144    0.6142    0.6565    0.5928    0.5546    0.5713    0.6247    0.6257    0.6126    0.6533];
            % b_param = [0.5899    0.4921    0.6374    0.5616    0.6330    0.5702    0.5963    0.6288    0.7326    0.4502    0.5968    0.6048    0.6458    0.5991    0.5918];
            % c_param = [42.6117   40.2665   37.8510   39.2992   37.6266   35.6608   36.3308   35.0988   35.2774   34.2218   37.7286   38.0913   38.2750   40.4708   44.0314];
            
            % Dodanie reguł TS w postaci liniowej
            for i = 1:15
                obj.linear_fis = addMF(obj.linear_fis, 'Y_fuzzy', 'linear', [a_param(i), b_param(i), c_param(i)]);
            end
            
            % Reguły Takagi-Sugeno: [inputMF, outputMF, weight]
            ruleList = [1 1 1 1 1;
                        1 2 2 1 1;
                        1 3 3 1 1;
                        2 1 4 1 1;
                        2 2 5 1 1;
                        2 3 6 1 1;
                        3 1 7 1 1;
                        3 2 8 1 1;
                        3 3 9 1 1;
                        4 1 10 1 1;
                        4 2 11 1 1;
                        4 3 12 1 1;
                        5 1 13 1 1;
                        5 2 14 1 1;
                        5 3 15 1 1];
            
            % Dodanie reguł do systemu
            obj.linear_fis = addRule(obj.linear_fis, ruleList);
        end

        function nonlinearFuzzy(obj)
            obj.nonlinear_fis.sigma_F1 = 25;
            obj.nonlinear_fis.sigma_FD = 15;
            obj.nonlinear_fis.rules_number = 6;
            obj.nonlinear_fis.F1_center = linspace(obj.F1_min, obj.F1_max, 3);
            obj.nonlinear_fis.FD_center = linspace(obj.FD_min, obj.FD_max, 2);

            obj.nonlinear_fis.a_param = [5.5084	5.2897	16.4147	18.9656	4.5280	7.4589];
            obj.nonlinear_fis.b_param = [0.6532	0.8794	0.6018	0.9663	2.3825	2.6390];
            obj.nonlinear_fis.c_param = [1.1386	8.8990	-11.4339	7.0581	3.7385	5.8944];

            % obj.nonlinear_fis.a_param = [3.1360    2.2899   10.7762   10.6849    4.6489    5.9102];
            % obj.nonlinear_fis.b_param = [0.3367    0.6262    1.2488    1.8362    1.1683    1.3966];
            % obj.nonlinear_fis.c_param = [24.0212   30.2372   29.2767   41.4704   39.3933   57.1661];
        end

        function testLinearModel(obj, U, a, b, obiekt, index)
            [Y_real, Y_lin] = obiekt.modifiedEuler(U, obiekt.kk);
            t = 0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp;
            
            Y_fuzzy = evalfis(obj.linear_fis, [U(1,:)', U(2,:)']);
            Y_out = zeros(1,obiekt.kk);
            
            for k = obiekt.delay+2:obiekt.kk
                if(U(1,k-(obiekt.delay+1)) ~= 0)
                    K = Y_fuzzy(k-(obiekt.delay+1)) / (U(1,k-(obiekt.delay+1)) + U(2, k-1));
                else
                    K = 1;
                end
                % Y_out(k) = -a * Y_out(k-1:-1:k-2)' + K * b * U(1, k-(obiekt.delay+1):-1:k-(obiekt.delay+2))' + K * b.F_D * U(2, k-1:-1:k-2)';
                Y_out(k) = -a * Y_out(k-1:-1:k-2)' + K * b * U(1, k-(obiekt.delay+1)) + K * b * U(2, k-1);
            end

            E_lin = sum((Y_real - Y_lin).^2) / obiekt.kk;
            E_out = sum((Y_real - Y_out).^2) / obiekt.kk;

            fprintf("\nHAMMERSTEIN LINEAR MODEL\n");
            fprintf("E_lin = %.3f\n", E_lin);
            fprintf("E_out = %.3f\n", E_out);
    
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

        function testNonlinearModel(obj, U, a, b, obiekt, index)
            gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));
            
            [Y_real, Y_lin] = obiekt.modifiedEuler(U, obiekt.kk);
            t = (0:length(U)-1) * obiekt.Tp; % Czas w sekundach (próbkowanie = 20s)
            
            Y_out = zeros(1,obiekt.kk);
            Y_fuzzy = zeros(1,obiekt.kk);
            
            for k = 1:obiekt.kk
                deg_u1 = [gaussmf_val(U(1,k), obj.nonlinear_fis.sigma_F1, obj.nonlinear_fis.F1_center(1)), ...
                    gaussmf_val(U(1,k), obj.nonlinear_fis.sigma_F1, obj.nonlinear_fis.F1_center(2)), ...
                    gaussmf_val(U(1,k), obj.nonlinear_fis.sigma_F1, obj.nonlinear_fis.F1_center(3))];
                deg_u2 = [gaussmf_val(U(2,k), obj.nonlinear_fis.sigma_FD, obj.nonlinear_fis.FD_center(1)), ...
                    gaussmf_val(U(2,k), obj.nonlinear_fis.sigma_FD, obj.nonlinear_fis.FD_center(2))];
                
                degrees_all(1) = deg_u1(1) * deg_u2(1);
                degrees_all(2) = deg_u1(1) * deg_u2(2);
                degrees_all(3) = deg_u1(2) * deg_u2(1);
                degrees_all(4) = deg_u1(2) * deg_u2(2);
                degrees_all(5) = deg_u1(3) * deg_u2(1);
                degrees_all(6) = deg_u1(3) * deg_u2(2);
                
                w = 0;
                output = 0;
                for i = 1:obj.nonlinear_fis.rules_number
                    output = output + degrees_all(i)*(obj.nonlinear_fis.a_param(i)*sinh(U(1,k)/22.5) + ...
                        obj.nonlinear_fis.b_param(i)*sinh(U(2,k)/7.5) + obj.nonlinear_fis.c_param(i));
                    w = w + degrees_all(i);
                end
                Y_fuzzy(k) = output / w;
            end

            for k = obiekt.delay+2:obiekt.kk
                if(U(1,k-(obiekt.delay+1)) ~= 0)
                    K = Y_fuzzy(k-(obiekt.delay+1)) / (U(1,k-(obiekt.delay+1)) + U(2, k-1));
                else
                    K = 1;
                end
                Y_out(k) = -a * Y_out(k-1:-1:k-2)' + K * b * U(1, k-(obiekt.delay+1)) + K * b * U(2, k-1);
            end

            E_lin = sum((Y_real - Y_lin).^2) / obiekt.kk;
            E_out = sum((Y_real - Y_out).^2) / obiekt.kk;
            fprintf("\nHAMMERSTEIN NONLINEAR MODEL\n");
            fprintf("E_lin = %.3f\n", E_lin);
            fprintf("E_out = %.3f\n", E_out);

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