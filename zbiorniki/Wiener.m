classdef Wiener < handle
    properties
        Y_min;
        Y_max;
        Y_threeCenter;
        Y_h
        h

        linear_fis
        nonlinear_fis
    end

    methods

        function obj = Wiener(a, b, obiekt)
            %Generacja danych sterujących i RK4

            U = [linspace(-45, 45, 100);
                linspace(-15, 15, 100)];

            [F1_grid, FD_grid] = meshgrid(U(1,:), U(2,:));  % kombinacje
            obj.h = ((F1_grid+obiekt.F_10 + FD_grid+obiekt.F_D0) ./ obiekt.alpha_2).^2 - obiekt.h_20;
            
            y = zeros(1, 100);
            for i = 1:100
                for j = 1:100
                    u = [ones(1,100) * U(1,j);
                        ones(1,100) * U(2,i)];
                    for k = 8:length(u)
                        y(k) = -a * y(k-1:-1:k-2)' + b * u(1, k-6:-1:k-7)' + b * u(2, k-1:-1:k-2)';
                    end
                    obj.Y_h(i,j) = y(end);
                end
            end

            obj.Y_min = min(min(obj.Y_h));
            obj.Y_max = max(max(obj.Y_h));
        end

        function linearFuzzy(obj)
            % % Dodanie reguł do systemu
            % obj.linear_fis = addRule(obj.linear_fis, ruleList);
            X = obj.Y_h(:);
            Y = obj.h(:);
            % Generowanie systemu rozmytego
            fis = genfis1([X Y], 5, 'gaussmf', 'linear');
            % Trening ANFIS
            options = anfisOptions('InitialFIS', fis, 'EpochNumber', 100, ...
                'DisplayANFISInformation', 0, 'DisplayErrorValues', 0);
            obj.linear_fis = anfis([X Y], options);
        end

        function nonlinearFuzzy(obj)
            X = sinh(obj.Y_h(:)/30);
            Y = obj.h(:);
            % Generowanie systemu rozmytego
            fis = genfis1([X Y], 3, 'gaussmf', 'linear');
            % Trening ANFIS
            options = anfisOptions('InitialFIS', fis, 'EpochNumber', 100, ...
                'DisplayANFISInformation', 0, 'DisplayErrorValues', 0);
            obj.nonlinear_fis = anfis([X Y], options);
        end

        function testLinearModel(obj, U, a, b, obiekt, index)
            [Y_real, Y_lin] = obiekt.modifiedEuler(U, obiekt.kk);
            t = 0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp;
           
            Y_out = zeros(1,obiekt.kk);
            Y_fuzzy = zeros(1,obiekt.kk);
            
            for k = obiekt.delay+3:obiekt.kk
                Y_out(k) = -a * Y_out(k-1:-1:k-2)' + b * U(1, k-(obiekt.delay+1):-1:k-(obiekt.delay+2))' + b * U(2, k-1:-1:k-2)';
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
            title('Porównanie wyjścia obiektu testowego i jego modeli - następniki liniowe');
            ylabel('y [cm]');
            xlabel('t [s]');
            grid on;

            % saveas(gcf, sprintf('D:/EiTI/MGR/raporty/raport_MGR/pictures/WienerLinearModel_%d.png', index));  % Zapisuje jako plik PNG
        end

        function testNonlinearModel(obj, U, a, b, obiekt, index)
            [Y_real, Y_lin] = obiekt.modifiedEuler(U, obiekt.kk);
            t = 0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp;
           
            Y_out = zeros(1,obiekt.kk);
            Y_fuzzy = zeros(1,obiekt.kk);
            
            for k = obiekt.delay+3:obiekt.kk
                Y_out(k) = -a * Y_out(k-1:-1:k-2)' + b * U(1, k-(obiekt.delay+1):-1:k-(obiekt.delay+2))' + b * U(2, k-1:-1:k-2)';
                Y_fuzzy(k) = evalfis(obj.nonlinear_fis, sinh(Y_out(k)/30));
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
            title('Porównanie wyjścia obiektu testowego i jego modeli - następniki nieliniowe');
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
