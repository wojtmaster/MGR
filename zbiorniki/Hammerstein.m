classdef Hammerstein < handle

    properties
        F1_min = -45;
        F1_max = 45;
        FD_min = -15;
        FD_max = 15;
        F1_grid
        FD_grid
        h

        linear_fis
        nonlinear_fis
    end

    methods
        function obj = Hammerstein(obiekt)
            U = [linspace(-45, 45, 100);
                linspace(-15, 15, 100)];
            [obj.F1_grid, obj.FD_grid] = meshgrid(U(1,:), U(2,:));  % kombinacje

            obj.h = ((obj.F1_grid+obiekt.F_10 + obj.FD_grid+obiekt.F_D0) ./ obiekt.alpha_2).^2 - obiekt.h_20;
        end

        function linearFuzzy(obj)
            X = [obj.F1_grid(:) obj.FD_grid(:)];
            Y = obj.h(:);
            % Generowanie systemu rozmytego
            fis = genfis1([X Y], 5, 'gaussmf', 'linear');
            % Trening ANFIS
            options = anfisOptions('InitialFIS', fis, 'EpochNumber', 100, ...
                'DisplayANFISInformation', 0, 'DisplayErrorValues', 0);
            obj.linear_fis = anfis([X Y], options);
        end

        function nonlinearFuzzy(obj)
            X = [sinh(obj.F1_grid(:)/45), sinh(obj.FD_grid(:)/15)];
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
            
            for k = obiekt.delay+3:obiekt.kk
                if(U(1,k-(obiekt.delay+1))+U(2, k-1) ~= 0)
                    Y_fuzzy =  evalfis(obj.linear_fis, [U(1,k-(obiekt.delay+1)), U(2,k-1)]);
                    K = Y_fuzzy / (U(1,k-(obiekt.delay+1)) + U(2, k-1));
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
    
            figure;
            plot(t, Y_real, 'b', t, Y_lin, 'g', t, Y_out, 'r', 'LineWidth', 1.5);

            legend({'model nieliniowy', ...
                    ['model liniowy' blanks(17) sprintf('E = %.3f', E_lin)], ...
                    ['model Hammersteina' blanks(5) sprintf('E = %.3f', E_out)]}, ...
                    'Location', 'best');

            title('Porównanie wyjścia obiektu testowego i jego modeli - następniki liniowe');
            ylabel('y [cm]');
            xlabel('t [s]');
            grid on;

            % saveas(gcf, sprintf('D:/EiTI/MGR/raporty/raport_MGR/pictures/HammersteinLinearModel_%d.png', index));  % Zapisuje jako plik PNG
        end

        function testNonlinearModel(obj, U, a, b, obiekt, index)
            [Y_real, Y_lin] = obiekt.modifiedEuler(U, obiekt.kk);
            t = 0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp; 
            
            Y_out = zeros(1,obiekt.kk);
            
            for k = obiekt.delay+3:obiekt.kk
                if(U(1,k-(obiekt.delay+1)) + U(2, k-1) ~= 0)
                    Y_fuzzy =  evalfis(obj.nonlinear_fis, [sinh(U(1,k-(obiekt.delay+1))/45), sinh(U(2,k-1)/15)]);
                    K = Y_fuzzy / (U(1,k-(obiekt.delay+1)) + U(2, k-1));
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
            title('Porównanie wyjścia obiektu testowego i jego modeli - następniki nieliniowe');
            ylabel('y [cm]');
            xlabel('t [s]');
            grid on;

            % saveas(gcf, sprintf('D:/EiTI/MGR/raporty/raport_MGR/pictures/HammersteinNonlinearModel_%d.png', index));  % Zapisuje jako plik PNG
        end

        function show_fuzzy(~, fis, s)
            figure;
            plotmf(fis, 'input', 1);
            title(sprintf('Funkcje przynależności u_1 - następniki %s', s));
            xlabel('u');
            ylabel('$\mu(u)$', 'Interpreter', 'latex');
            grid on;

            figure;
            plotmf(fis, 'input', 2);
            title(sprintf('Funkcje przynależności u_2 - następniki %s', s));
            xlabel('u');
            ylabel('$\mu(u)$', 'Interpreter', 'latex');
            grid on;

            % saveas(gcf, sprintf('D:/EiTI/MGR/raporty/raport_MGR/pictures/HammersteinfuzzySets_%s.png', s));  % Zapisuje jako plik PNG
        end
    end
end