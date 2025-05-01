classdef Hammerstein < handle

    properties
        U_min = -45;
        U_max = 45;
        U_threeCenter;

        linear_fis
        nonlinear_fis
    end

    methods
        function obj = Hammerstein()
        end

        function linearFuzzy(obj)
            U_center = linspace(obj.U_min, obj.U_max, 5);
            obj.linear_fis = sugfis('Name', 'Hammerstein', 'Type', 'sugeno');
            obj.linear_fis = addInput(obj.linear_fis, [obj.U_min obj.U_max], 'Name', 'U');
            
            % Definiowanie funkcji przynależności (gaussmf)
            for i = 1:length(U_center) 
                obj.linear_fis = addMF(obj.linear_fis, 'U', 'gaussmf', [12, U_center(i)]);
            end
            
            % Definiowanie wyjścia i początkowych następników (a_i * u + b_i)
            obj.linear_fis = addOutput(obj.linear_fis, [obj.U_min obj.U_max], 'Name', 'U_fuzzy');
            
            % Optymalne parametry a: 
            a_param = [0.3647    0.4227    0.4844    0.5657    0.6796];
            % Optymalne parametry b: 
            b_param = [30.4365   33.2982   36.0644   38.2423   38.0915];
            
            % Dodanie reguł TS w postaci liniowej
            for i = 1:length(U_center)
                obj.linear_fis = addMF(obj.linear_fis, 'U_fuzzy', 'linear', [a_param(i), b_param(i)]);
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
            obj.U_threeCenter = linspace(obj.U_min, obj.U_max, 3);
            obj.nonlinear_fis.sigma = 25;
            % Optymalne parametry a: 
            obj.nonlinear_fis.a_param = [0.4039    1.3416    0.3739];
            % Optymalne parametry b: 
            obj.nonlinear_fis.b_param = [17.1217   22.6325   11.9320];
            % Optymalne parametry c: 
            obj.nonlinear_fis.c_param = [13.8112   34.3922   65.9774];
        end

        function testLinearModel(obj, U, a, b, obiekt, index)
            [Y_real, Y_lin] = obiekt.rk4(U, obiekt.kk);
            t = 0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp;
            
            Y_fuzzy = evalfis(obj.linear_fis, U(1,:));
            Y_0 = evalfis(obj.linear_fis, -10);
            Y_out = zeros(1,obiekt.kk);
            
            for k = 8:obiekt.kk
                if(U(1,k-(obiekt.delay+2)) ~= 0)
                    K = (Y_fuzzy(k-(obiekt.delay+2)) - Y_0) / U(1,k-(obiekt.delay+2));
                else
                    K = 1;
                end
                Y_out(k) = -a * Y_out(k-1:-1:k-2)' + K * b * U(1, k-(obiekt.delay+1):-1:k-(obiekt.delay+2))' + K * b * U(2, k-1:-1:k-2)';
            end

            E_lin = sum((Y_real - Y_lin).^2) / obiekt.kk;
            E_out = sum((Y_real - Y_out).^2) / obiekt.kk;

            fprintf("\nHAMMERSTEIN LINEAR MODEL\n");
            fprintf("E_lin = %.3f\n", E_lin);
            fprintf("E_out = %.3f\n", E_out);
    
            % figure;
            % plot(t, Y_real, 'b', t, Y_lin, 'g', t, Y_out, 'r', 'LineWidth', 1.5);
            % 
            % legend({'model nieliniowy', ...
            %         ['model liniowy' blanks(17) sprintf('E = %.3f', E_lin)], ...
            %         ['model Hammersteina' blanks(5) sprintf('E = %.3f', E_out)]}, ...
            %         'Location', 'best');
            % 
            % title('Porównanie wyjścia obiektu testowego i jego modeli');
            % ylabel('y [cm]');
            % xlabel('t [s]');
            % grid on;

            % saveas(gcf, sprintf('D:/EiTI/MGR/raporty/raport_MGR/pictures/HammersteinLinearModel_%d.png', index));  % Zapisuje jako plik PNG
        end

        function testNonlinearModel(obj, U, a, b, obiekt, index)
            gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));

            [Y_real, Y_lin] = obiekt.rk4(U, obiekt.kk);
            t = 0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp;
            
            Y_out = zeros(1,obiekt.kk);
            Y_fuzzy = zeros(1,obiekt.kk);
            
            output = 0;
            w = 0;
            degrees = [gaussmf_val(0, obj.nonlinear_fis.sigma, obj.U_threeCenter(1)), ...
                       gaussmf_val(0, obj.nonlinear_fis.sigma, obj.U_threeCenter(2)), ...
                       gaussmf_val(0, obj.nonlinear_fis.sigma, obj.U_threeCenter(3))];
            for i = 1:length(degrees)
                output = output +  degrees(i) * (obj.nonlinear_fis.a_param(i)*sinh(0/obj.nonlinear_fis.b_param(i)) + obj.nonlinear_fis.c_param(i));
                w = w + degrees(i);
            end
            Y_0 = output / w;
            
            for k = 1:7
                output = 0;
                w = 0;
                degrees = [gaussmf_val(U(1,k), obj.nonlinear_fis.sigma, obj.U_threeCenter(1)), ...
                           gaussmf_val(U(1,k), obj.nonlinear_fis.sigma, obj.U_threeCenter(2)), ...
                           gaussmf_val(U(1,k), obj.nonlinear_fis.sigma, obj.U_threeCenter(3))];
                for i = 1:length(degrees)
                    output = output +  degrees(i) * (obj.nonlinear_fis.a_param(i)*sinh(U(1,k)/obj.nonlinear_fis.b_param(i)) + obj.nonlinear_fis.c_param(i));
                    w = w + degrees(i);
                end
                Y_fuzzy(k) = output / w;
            end
            
            for k = 8:obiekt.kk
                output = 0;
                w = 0;
                degrees = [gaussmf_val(U(1,k), obj.nonlinear_fis.sigma, obj.U_threeCenter(1)), ...
                           gaussmf_val(U(1,k), obj.nonlinear_fis.sigma, obj.U_threeCenter(2)), ...
                           gaussmf_val(U(1,k), obj.nonlinear_fis.sigma, obj.U_threeCenter(3))];
                for i = 1:length(degrees)
                    output = output +  degrees(i) * (obj.nonlinear_fis.a_param(i)*sinh(U(1,k)/obj.nonlinear_fis.b_param(i)) + obj.nonlinear_fis.c_param(i));
                    w = w + degrees(i);
                end
                Y_fuzzy(k) = output / w;
            
                if(U(1,k-(obiekt.delay+2)) ~= 0)
                    K = (Y_fuzzy(k-(obiekt.delay+2)) - Y_0) / U(1,k-(obiekt.delay+2));
                else
                    K = 1;
                end
                Y_out(k) = -a * Y_out(k-1:-1:k-2)' + K * b * U(1, k-(obiekt.delay+1):-1:k-(obiekt.delay+2))' + K * b * U(2, k-1:-1:k-2)';
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