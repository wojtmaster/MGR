classdef Wiener < handle
    properties
        Y_min;
        Y_max;
        Y_threeCenter;

        linear_fis
        nonlinear_fis
    end

    methods

        function obj = Wiener(a, b)
            %Generacja danych sterujących i RK4
            U = [linspace(-45, 45, 100);
                linspace(-15, 15, 100)];
            
            h = zeros(100, 100);
            y = zeros(1, 100);
            for i = 1:100
                for j = 1:100
                    u = [ones(1,100) * U(1,j);
                        ones(1,100) * U(2,i)];
                    for k = 8:length(u)
                        y(k) = -a * y(k-1:-1:k-2)' + b * u(1, k-6) + b * u(2, k-1);
                    end
                    h(i,j) = y(end);
                end
            end

            obj.Y_min = min(min(h));
            obj.Y_max = max(max(h));
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
            
            a_param = [	0.5079	0.5508	0.5985	0.6438	0.6801];
            b_param = [	3.5118	0.5148	-0.1804	0.6402	4.2028];

            % a_param = [0.5409    0.5918    0.6604    0.7429    0.8309];
            % b_param = [41.5737   37.8453   36.0261   34.0338   31.3233];
            
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
            obj.nonlinear_fis.Y_center = linspace(obj.Y_min, obj.Y_max, 3);
            obj.nonlinear_fis.sigma = 30;
            obj.nonlinear_fis.rules_number = 3;
            obj.nonlinear_fis.sigma = 30;

            obj.nonlinear_fis.a_param = [5.7627	16.2006	12.0295];
            obj.nonlinear_fis.b_param = [-10.9354 -0.5527 15.2669];
            
            % obj.nonlinear_fis.a_param = [3.6419   13.1436    7.7333];
            % obj.nonlinear_fis.b_param = [24.4050   35.2421   53.3554];
        end

        function testLinearModel(obj, U, a, b, obiekt, index)
            [Y_real, Y_lin] = obiekt.modifiedEuler(U, obiekt.kk);
            t = 0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp;
           
            Y_out = zeros(1,obiekt.kk);
            Y_fuzzy = zeros(1,obiekt.kk);
            
            for k = obiekt.delay+2:obiekt.kk
                Y_out(k) = -a * Y_out(k-1:-1:k-2)' + b * U(1, k-(obiekt.delay+1)) + b * U(2, k-1);
                Y_fuzzy(k) = evalfis(obj.linear_fis, Y_out(k));
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

            [Y_real, Y_lin] = obiekt.modifiedEuler(U, obiekt.kk);
            t = (0:length(U)-1) * obiekt.Tp; % Czas w sekundach (próbkowanie = 20s)
            
            Y_out = zeros(1, obiekt.kk);
            Y_fuzzy = zeros(1, obiekt.kk);
            
            for k = obiekt.delay+2:obiekt.kk
                Y_out(k) = -a * Y_out(k-1:-1:k-2)' + b * U(1, k-(obiekt.delay+1)) + b * U(2, k-1);
                w = 0;
                output = 0;
                degrees = [gaussmf_val(Y_out(k), obj.nonlinear_fis.sigma, obj.nonlinear_fis.Y_center(1)), ...
                           gaussmf_val(Y_out(k), obj.nonlinear_fis.sigma, obj.nonlinear_fis.Y_center(2)), ...
                           gaussmf_val(Y_out(k), obj.nonlinear_fis.sigma, obj.nonlinear_fis.Y_center(3))];
                for i = 1:obj.nonlinear_fis.rules_number
                    output = output +  degrees(i) * (obj.nonlinear_fis.a_param(i)*sinh(Y_out(k)/36) + obj.nonlinear_fis.b_param(i));
                    w = w + degrees(i);
                end
                Y_fuzzy(k) = output / w;
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
