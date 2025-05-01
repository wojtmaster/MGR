classdef Wiener < handle
    properties
        Y_min;
        Y_max;
        Y_threeCenter;

        linear_fis
        nonlinear_fis
    end

    methods

        function obj = Wiener(obiekt, a, b)
            %Generacja danych sterujących i RK4
            U = linspace(-45, 45, obiekt.kk);
            
            h = zeros(1, obiekt.kk);
            y = zeros(1, 100);
            for i = 1:obiekt.kk
                u = [ones(1,100) * U(i);
                    zeros(1,100)];
                for k = 8:length(u)
                    y(k) = -a * y(k-1:-1:k-2)' + b * u(1, k-(obiekt.delay+1):-1:k-(obiekt.delay+2))' + b * u(2, k-1:-1:k-2)';
                end
                h(i) = y(end);
            end

            obj.Y_min = min(h);
            obj.Y_max = max(h);
        end

        function linearFuzzy(obj)
            Y_center = linspace(obj.Y_min, obj.Y_max, 5); % Środek zbiorów
            
            obj.linear_fis = sugfis('Name', 'Wiener', 'Type', 'sugeno');
            obj.linear_fis = addInput(obj.linear_fis, [obj.Y_min obj.Y_max], 'Name', 'Y');
            
            % Definiowanie funkcji przynależności (gaussmf)
            for i = 1:length(Y_center) 
                obj.linear_fis = addMF(obj.linear_fis, 'Y', 'gaussmf', [12, Y_center(i)]);
            end
            
            % Definiowanie wyjścia i początkowych następników (a_i * u + b_i)
            obj.linear_fis = addOutput(obj.linear_fis, [obj.Y_min obj.Y_max], 'Name', 'Y_fuzzy');
            
            % Optymalne parametry a: 
            a_param = [0.3647    0.4227    0.4844    0.5657    0.6796];
            % Optymalne parametry b: 
            b_param = [30.4365   33.2982   36.0644   38.2423   38.0915];
            
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
            obj.Y_threeCenter = linspace(obj.Y_min, obj.Y_max, 3);
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
           
            Y_out= zeros(1,obiekt.kk);
            Y_fuzzy = zeros(1,obiekt.kk);
            
            for k = 8:obiekt.kk
                Y_out(k) = -a * Y_out(k-1:-1:k-2)' + b * U(1, k-(obiekt.delay+1):-1:k-(obiekt.delay+2))' + b * U(2, k-1:-1:k-2)';
                Y_fuzzy(k) = evalfis(obj.linear_fis, Y_out(k)) - obiekt.h_20;
            end

            E_lin = sum((Y_real - Y_lin).^2) / obiekt.kk;
            E_out = sum((Y_real - Y_fuzzy).^2) / obiekt.kk;
            fprintf("\nWIENER LINEAR MODEL\n");
            fprintf("E_lin = %.3f\n", E_lin);
            fprintf("E_out = %.3f\n", E_out);

            figure;
            plot(t, Y_real, 'b', t, Y_lin, 'g', t, Y_fuzzy, 'r', 'LineWidth', 1.5);
            legend('model nieliniowy', ...
                sprintf('model liniowy \t E = %.3f', E_lin), ...
                sprintf('model Wiener \t E = %.3f', E_out), 'Location', 'best');
            title('Porównanie wyjścia obiektu testowego i jego modeli');
            ylabel('y [cm]');
            xlabel('t [s]');
            grid on;

            % saveas(gcf, sprintf('D:/EiTI/MGR/raporty/raport_MGR/pictures/WienerLinearModel_%d.png', index));  % Zapisuje jako plik PNG
        end

        function testNonlinearModel(obj, U, a, b, obiekt, index)
            gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));

            [Y_real, Y_lin] = obiekt.rk4(U, obiekt.kk);
            t = 0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp;

            Y_out = zeros(1,obiekt.kk);
            Y_fuzzy = zeros(1,obiekt.kk);
            
            for k = 8:obiekt.kk
                Y_out(k) = -a * Y_out(k-1:-1:k-2)' + b * U(1, k-(obiekt.delay+1):-1:k-(obiekt.delay+2))' + b * U(2, k-1:-1:k-2)';
                output = 0;
                w = 0;
                degrees = [gaussmf_val(Y_out(k), obj.nonlinear_fis.sigma, obj.Y_threeCenter(1)), ...
                           gaussmf_val(Y_out(k), obj.nonlinear_fis.sigma, obj.Y_threeCenter(2)), ...
                           gaussmf_val(Y_out(k), obj.nonlinear_fis.sigma, obj.Y_threeCenter(3))];
                for i = 1:length(degrees)
                    output = output +  degrees(i) * (obj.nonlinear_fis.a_param(i)*sinh(Y_out(k)/obj.nonlinear_fis.b_param(i)) + obj.nonlinear_fis.c_param(i));
                    w = w + degrees(i);
                end
                Y_fuzzy(k) = output / w - obiekt.h_20;
            end
            
            E_lin = sum((Y_real - Y_lin).^2) / obiekt.kk;
            E_out = sum((Y_real - Y_fuzzy).^2) / obiekt.kk;
            fprintf("\nWIENER NONLINEAR MODEL\n");
            fprintf("E_lin = %.3f\n", E_lin);
            fprintf("E_out = %.3f\n", E_out);

            figure;
            plot(t, Y_real, 'b', t, Y_lin, 'g', t, Y_fuzzy, 'r', 'LineWidth', 1.5);
            legend('model nieliniowy', ...
                sprintf('model liniowy \t E = %.3f', E_lin), ...
                sprintf('model Wiener \t E = %.3f', E_out), 'Location', 'best');
            title('Porównanie wyjścia obiektu testowego i jego modeli');
            ylabel('y [cm]');
            xlabel('t [s]');
            grid on;

            % saveas(gcf, sprintf('D:/EiTI/MGR/raporty/raport_MGR/pictures/WienerNonlinearModel_%d.png', index));  % Zapisuje jako plik PNG
        end
        
        function show_fuzzy(~, fis, s)
            figure;
            plotmf(fis, 'input', 1);
            title(sprintf('Funkcje przynależności dla Y - następniki %s', s));
            xlabel('y');
            ylabel('$\mu(y)$', 'Interpreter', 'latex');
            grid on;

            % saveas(gcf, sprintf('D:/EiTI/MGR/raporty/raport_MGR/pictures/WienerfuzzySets_%s.png', s));  % Zapisuje jako plik PNG
        end

    end
end
