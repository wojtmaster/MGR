classdef Wiener < handle
    properties
        Y_h
        Y_pH
        h
        pH
        Y_min = -30;
        Y_max = 30;
        U_min = -15;
        U_max= 15;

        linear_fis
        nonlinear_fis
    end

    methods

        function obj = Wiener(a_h, a_pH, b_h, b_pH, obiekt)
            U = [linspace(1.6, 31.6, 100)
                linspace(0.6, 30.6, 100)];
            [Q1_grid, Q3_grid] = meshgrid(U(1,:), U(2,:));
            obj.h = ((Q1_grid + obiekt.Q_20 + Q3_grid) / obiekt.C_V).^2 - obiekt.h_0;
            Wa4 = (obiekt.W_a1*Q1_grid + obiekt.W_a2*obiekt.Q_20 + obiekt.W_a3*Q3_grid)./(Q1_grid+obiekt.Q_20+Q3_grid);
            Wb4 = (obiekt.W_b1*Q1_grid + obiekt.W_b2*obiekt.Q_20 + obiekt.W_b3*Q3_grid)./(Q1_grid+obiekt.Q_20+Q3_grid);
            
            for i = 1:100
                for j = 1:100
                    obj.pH(i,j) = obiekt.pH_calc(Wa4(i,j), Wb4(i,j)) - obiekt.pH_0;
                end
            end

            U = [linspace(obj.U_min, obj.U_max, 100);
                linspace(obj.U_min, obj.U_max, 100)];

            [Q1_grid, Q3_grid] = meshgrid(U(1,:), U(2,:));

            obj.Y_h = zeros(100, 100);
            obj.Y_pH = zeros(100, 100);
            
            y_h = zeros(1, 100);
            y_1 = zeros(1, 100);
            y_2 = zeros(1, 100);
            y_pH = zeros(1, 100);
            
            for i = 1:100
                for j = 1:100
                    u = [ones(1,100) * U(1,j);
                        ones(1,100) * U(2,i)];
                    for k = 2:length(u)
                        y_h(k) = - a_h*y_h(k-1) + b_h*u(1, k-1) + b_h*u(2, k-1);
            
                        if k >= 4
                            y_1(k) = -a_pH.Q1 * y_1(k-1:-1:k-2)' + b_pH.Q1*u(1,k-1:-1:k-2)';
                            y_2(k) = -a_pH.Q3 * y_2(k-1:-1:k-2)' + b_pH.Q3*u(2,k-1:-1:k-2)';
                            y_pH(k) = y_1(k) + y_2(k);
                        end
                    end
                    obj.Y_h(i,j) = y_h(end);
                    obj.Y_pH(i,j) = y_pH(end);
                end
            end

            figure;
            surf(Q1_grid, Q3_grid, obj.Y_h);

            figure;
            surf(Q1_grid, Q3_grid, obj.Y_pH);
        end

        function linearFuzzy(obj)
            % -------------------- h --------------------%
            X = obj.Y_h';
            X = X(:);
            Y = obj.h';
            Y = Y(:);
            % Generowanie systemu rozmytego
            fis = genfis1([X Y], 5, 'gaussmf', 'linear');
            % Trening ANFIS
            options = anfisOptions('InitialFIS', fis, 'EpochNumber', 100, ...
                'DisplayANFISInformation', 0, 'DisplayErrorValues', 0);
            obj.linear_fis.h = anfis([X Y], options);

            X = obj.Y_pH';
            X = X(:);
            Y = obj.pH';
            Y = Y(:);
            % Generowanie systemu rozmytego
            fis = genfis1([X Y], 12, 'gaussmf', 'linear');
            % Trening ANFIS
            options = anfisOptions('InitialFIS', fis, 'EpochNumber', 200, ...
                'DisplayANFISInformation', 0, 'DisplayErrorValues', 0);
            obj.linear_fis.pH = anfis([X Y], options);
        end

        function nonlinearFuzzy(obj)
            % ------------------ h ------------------ %
            X = obj.Y_h';
            X = sinh(X(:)/30);
            Y = obj.h;
            Y = Y(:);
            % Generowanie systemu rozmytego
            fis = genfis1([X Y], 3, 'gaussmf', 'linear');
            % Trening ANFIS
            options = anfisOptions('InitialFIS', fis, 'EpochNumber', 100, ...
                'DisplayANFISInformation', 0, 'DisplayErrorValues', 0);
            obj.nonlinear_fis.h = anfis([X Y], options);

            X = obj.Y_pH';
            X = tanh(X(:)/15);
            Y = obj.pH';
            Y = Y(:);
            % Generowanie systemu rozmytego
            fis = genfis1([X Y], 4, 'gaussmf', 'linear');
            % Trening ANFIS
            options = anfisOptions('InitialFIS', fis, 'EpochNumber', 200, ...
                'DisplayANFISInformation', 0, 'DisplayErrorValues', 0);
            obj.nonlinear_fis.pH = anfis([X Y], options);
        end

        function testLinearModel(obj, U, a_h, a_pH, b_h, b_pH, obiekt, type, index)
            [Y_real, Y_lin] = obiekt.modifiedEuler(U, obiekt.kk); % Symulacja rzeczywistego układu
            t = 0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp; 

            v = zeros(1, obiekt.kk);
            Y_out = zeros(1, obiekt.kk);
            Y_fuzzy = zeros(1, obiekt.kk);

            % -------------------- LINEAR H -------------------- %
            if(strcmp(type, 'h'))
                for k = 2:obiekt.kk
                    Y_out(k) = - a_h*Y_out(k-1) + b_h*U(1,k-1) + b_h * U(2, k-1) + b_h * U(3,k-1);
                    Y_fuzzy(k) = evalfis(obj.linear_fis.h, Y_out(k));
                end
                E_lin = sum((Y_real(1, :) - Y_lin(1, :)).^2) / obiekt.kk;
                E_out = sum((Y_real(1, :) - Y_fuzzy).^2) / obiekt.kk;
                fprintf("\nWIENER LINEAR MODEL\n");
                fprintf("E_lin = %.3f\n", E_lin);
                fprintf("E_out = %.3f\n", E_out);
    
                figure;
                plot(t, Y_real(1, :), 'b', t, Y_lin(1, :), 'g', t, Y_fuzzy, 'r', 'LineWidth', 1.5);
                
                legend({'model nieliniowy', ...
                        ['model liniowy' blanks(17) sprintf('E = %.3f', E_lin)], ...
                        ['model Wiener' blanks(5) sprintf('E = %.3f', E_out)]}, ...
                        'Location', 'best');
    
                title('Porównanie wyjścia obiektu testowego i jego modeli');
                ylabel('h [cm]');
                xlabel('t [s]');
                grid on;

            % -------------------- LINEAR PH -------------------- %
            else
                y_1 = zeros(1,obiekt.kk);
                y_2 = zeros(1,obiekt.kk);
                y_3 = zeros(1,obiekt.kk);
                for k = 4:obiekt.kk
                    y_1(k) = -a_pH.Q1 * y_1(k-1:-1:k-2)' + b_pH.Q1*U(1,k-1:-1:k-2)';
                    y_2(k) = -a_pH.Q2 * y_2(k-1:-1:k-2)' + b_pH.Q2*U(2,k-1:-1:k-2)';
                    y_3(k) = -a_pH.Q3 * y_3(k-1:-1:k-2)' + b_pH.Q3*U(3,k-1:-1:k-2)';
                
                    Y_out(k) = y_1(k) + y_2(k) + y_3(k);
                    
                    Y_fuzzy(k) = evalfis(obj.linear_fis.pH, Y_out(k));
                end
                E_lin = sum((Y_real(2, :) - Y_lin(2, :)).^2) / obiekt.kk;
                E_out = sum((Y_real(2, :) - Y_fuzzy).^2) / obiekt.kk;
                fprintf("\nWIENER LINEAR MODEL\n");
                fprintf("E_lin = %.3f\n", E_lin);
                fprintf("E_out = %.3f\n", E_out);
    
                figure;
                plot(t, Y_real(2, :), 'b', t, Y_lin(2, :), 'g', t, Y_fuzzy, 'r', 'LineWidth', 1.5);
                
                legend({'model nieliniowy', ...
                        ['model liniowy' blanks(17) sprintf('E = %.3f', E_lin)], ...
                        ['model Wiener' blanks(5) sprintf('E = %.3f', E_out)]}, ...
                        'Location', 'best');
    
                title('Porównanie wyjścia obiektu testowego i jego modeli');
                ylabel('pH');
                xlabel('t [s]');
                grid on;
            end

            % saveas(gcf, sprintf('D:/EiTI/MGR/raporty/raport_MGR/pictures/WienerLinearModel_%d.png', index));  % Zapisuje jako plik PNG
        end

        function testNonlinearModel(obj, U, a_h, a_pH, b_h, b_pH, obiekt, type, index)
            [Y_real, Y_lin] = obiekt.modifiedEuler(U, obiekt.kk); % Symulacja rzeczywistego układu
            t = 0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp;

            Y_out = zeros(1,obiekt.kk);
            Y_fuzzy = zeros(1,obiekt.kk);

            % -------------------- NONLINEAR H -------------------- %
            if(strcmp(type, 'h'))
                for k = 2:obiekt.kk
                    Y_out(k) = - a_h*Y_out(k-1) + b_h*U(1,k-1) + b_h * U(2, k-1) + b_h * U(3,k-1);
                    Y_fuzzy(k) = evalfis(obj.nonlinear_fis.h, sinh(Y_out(k)/30));
                end
                E_lin = sum((Y_real(1, :) - Y_lin(1, :)).^2) / obiekt.kk;
                E_out = sum((Y_real(1, :) - Y_fuzzy).^2) / obiekt.kk;
                fprintf("\nWIENER LINEAR MODEL\n");
                fprintf("E_lin = %.3f\n", E_lin);
                fprintf("E_out = %.3f\n", E_out);
    
                figure;
                plot(t, Y_real(1, :), 'b', t, Y_lin(1, :), 'g', t, Y_fuzzy, 'r', 'LineWidth', 1.5);
                
                legend({'model nieliniowy', ...
                        ['model liniowy' blanks(17) sprintf('E = %.3f', E_lin)], ...
                        ['model Wiener' blanks(5) sprintf('E = %.3f', E_out)]}, ...
                        'Location', 'best');
    
                title('Porównanie wyjścia obiektu testowego i jego modeli');
                ylabel('h [cm]');
                xlabel('t [s]');
                grid on;

            % -------------------- NONLINEAR PH -------------------- %
            else
                y_1 = zeros(1,obiekt.kk);
                y_2 = zeros(1,obiekt.kk);
                y_3 = zeros(1,obiekt.kk);
                for k = 4:obiekt.kk
                    y_1(k) = -a_pH.Q1 * y_1(k-1:-1:k-2)' + b_pH.Q1*U(1,k-1:-1:k-2)';
                    y_2(k) = -a_pH.Q2 * y_2(k-1:-1:k-2)' + b_pH.Q2*U(2,k-1:-1:k-2)';
                    y_3(k) = -a_pH.Q3 * y_3(k-1:-1:k-2)' + b_pH.Q3*U(3,k-1:-1:k-2)';
                
                    Y_out(k) = y_1(k) + y_2(k) + y_3(k);
                    Y_fuzzy(k) = evalfis(obj.nonlinear_fis.pH, tanh(Y_out(k)/15));
                end
                E_lin = sum((Y_real(2, :) - Y_lin(2, :)).^2) / obiekt.kk;
                E_out = sum((Y_real(2, :) - Y_fuzzy).^2) / obiekt.kk;
                fprintf("\nWIENER NONLINEAR MODEL\n");
                fprintf("E_lin = %.3f\n", E_lin);
                fprintf("E_out = %.3f\n", E_out);
    
                figure;
                plot(t, Y_real(2, :), 'b', t, Y_lin(2, :), 'g', t, Y_fuzzy, 'r', 'LineWidth', 1.5);
                
                legend({'model nieliniowy', ...
                        ['model liniowy' blanks(17) sprintf('E = %.3f', E_lin)], ...
                        ['model Wiener' blanks(5) sprintf('E = %.3f', E_out)]}, ...
                        'Location', 'best');
    
                title('Porównanie wyjścia obiektu testowego i jego modeli');
                ylabel('pH');
                xlabel('t [s]');
                grid on;
            end

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
