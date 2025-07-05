classdef Hammerstein < handle

    properties
        U_min = -15;
        U_max = 15;
        Q1_grid
        Q3_grid
        h
        Wa4
        Wb4
        pH

        linear_fis
        nonlinear_fis
    end

    methods
        function obj = Hammerstein(obiekt)
            U = [linspace(1.6, 31.6, 100)
                linspace(0.6, 30.6, 100)];
            [obj.Q1_grid, obj.Q3_grid] = meshgrid(U(1,:), U(2,:));
            obj.h = ((obj.Q1_grid + obiekt.Q_20 + obj.Q3_grid) / obiekt.C_V).^2 - obiekt.h_0;
            obj.Wa4 = (obiekt.W_a1*obj.Q1_grid + obiekt.W_a2*obiekt.Q_20 + obiekt.W_a3*obj.Q3_grid)./(obj.Q1_grid+obiekt.Q_20+obj.Q3_grid);
            obj.Wb4 = (obiekt.W_b1*obj.Q1_grid + obiekt.W_b2*obiekt.Q_20 + obiekt.W_b3*obj.Q3_grid)./(obj.Q1_grid+obiekt.Q_20+obj.Q3_grid);
            
            for i = 1:100
                for j = 1:100
                    obj.pH(i,j) = obiekt.pH_calc(obj.Wa4(i,j), obj.Wb4(i,j)) - obiekt.pH_0;
                end
            end

            obj.Wa4 = obj.Wa4 - obiekt.W_a40;
            obj.Wb4 = obj.Wb4 - obiekt.W_b40;

            U = [linspace(obj.U_min, obj.U_max, 100);
                linspace(obj.U_min, obj.U_max, 100)];
            [obj.Q1_grid, obj.Q3_grid] = meshgrid(U(1,:), U(2,:));  % kombinacje
        end

        function linearFuzzy(obj)
            % --------------- LINEAR H --------------- %
            X1 = obj.Q1_grid';  X2 = obj.Q3_grid';
            X = [X1(:) X2(:)];
            
            Y = obj.h';
            Y = Y(:);

            % Generowanie systemu rozmytego
            fis = genfis1([X Y], 5, 'gaussmf', 'linear');
            % Trening ANFIS
            options = anfisOptions('InitialFIS', fis, 'EpochNumber', 100, ...
                'DisplayANFISInformation', 0, 'DisplayErrorValues', 0);
            obj.linear_fis.h = anfis([X Y], options);

            % --------------- LINEAR PH --------------- %
            X = [X1(:) X2(:)];
           
            Y_Wa4 = obj.Wa4';
            Y_Wa4 = Y_Wa4(:);
            
            Y_Wb4 = obj.Wb4';
            Y_Wb4 = Y_Wb4(:);

            % Generowanie systemu rozmytego
            fis_Wa4 = genfis1([X Y_Wa4], 5, 'gaussmf', 'linear');
            fis_Wb4 = genfis1([X Y_Wb4], 5, 'gaussmf', 'linear');
            
            % Trening ANFIS
            options_1 = anfisOptions('InitialFIS', fis_Wa4, 'EpochNumber', 100, ...
                'DisplayANFISInformation', 0, 'DisplayErrorValues', 0);
            options_2 = anfisOptions('InitialFIS', fis_Wb4, 'EpochNumber', 100, ...
                'DisplayANFISInformation', 0, 'DisplayErrorValues', 0);
            obj.linear_fis.Wa4 = anfis([X Y_Wa4], options_1);
            obj.linear_fis.Wb4 = anfis([X Y_Wb4], options_2);
        end

        function nonlinearFuzzy(obj)
            % --------------- NONLINEAR H --------------- %
            X1 = obj.Q1_grid';  X2 = obj.Q3_grid';
            X = [sinh(X1(:)/15) sinh(X2(:)/15)];

            Y = obj.h';
            Y = Y(:);
            
            % Generowanie systemu rozmytego
            fis = genfis1([X Y], 3, 'gaussmf', 'linear');
            % Trening ANFIS
            options = anfisOptions('InitialFIS', fis, 'EpochNumber', 100, ...
                'DisplayANFISInformation', 0, 'DisplayErrorValues', 0);
            obj.nonlinear_fis.h = anfis([X Y], options);
        
            % % --------------- NONLINEAR PH --------------- %
            % Dane wejściowe
            X_Wa4 = [tanh(X1(:)/15) tanh(X2(:)/15)];
            X_Wb4 = [sinh(X1(:)/15) sinh(X2(:)/15)];
            
            Y_Wa4 = obj.Wa4';
            Y_Wa4 = Y_Wa4(:);
            
            Y_Wb4 = obj.Wb4';
            Y_Wb4 = Y_Wb4(:);
            
            % Generowanie systemu rozmytego
            fis_Wa4 = genfis1([X_Wa4 Y_Wa4], 3, 'gaussmf', 'linear');
            fis_Wb4 = genfis1([X_Wb4 Y_Wb4], 3, 'gaussmf', 'linear');
            
            % Trening ANFIS
            options_1 = anfisOptions('InitialFIS', fis_Wa4, 'EpochNumber', 100, ...
                'DisplayANFISInformation', 0, 'DisplayErrorValues', 0);
            options_2 = anfisOptions('InitialFIS', fis_Wb4, 'EpochNumber', 100, ...
                'DisplayANFISInformation', 0, 'DisplayErrorValues', 0);
            obj.nonlinear_fis.Wa4 = anfis([X_Wa4 Y_Wa4], options_1);
            obj.nonlinear_fis.Wb4 = anfis([X_Wb4 Y_Wb4], options_2);
        end

        function testLinearModel(obj, U, a_h, a_pH, b_h, b_pH, obiekt, type, index)
            [Y_real, Y_lin] = obiekt.modifiedEuler(U, obiekt.kk);
            t = 0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp; 
            Y_out = zeros(1, obiekt.kk);
            
            % --------------- LINEAR H --------------- %
            if(strcmp(type, 'h'))
                Y_fuzzy = evalfis(obj.linear_fis.h, [U(1,:)', U(3,:)']);
                for k = 2:obiekt.kk
                    % Y_fuzzy(k) = ((obiekt.Q_10+obiekt.Q_20+obiekt.Q_30+U(1,k-1)+U(3,k-1)) / obiekt.C_V).^2 - obiekt.h_0;
                    Y_out(k) = - a_h*Y_out(k-1) + b_h*Y_fuzzy(k);
                end

                E_lin = sum((Y_real(1, :) - Y_lin(1, :)).^2) / obiekt.kk;
                E_out = sum((Y_real(1, :) - Y_out).^2) / obiekt.kk;
                fprintf("\nHAMMERSTEIN LINEAR MODEL\n");
                fprintf("E_lin = %.3f\n", E_lin);
                fprintf("E_out = %.3f\n", E_out);
    
                figure;
                plot(t, Y_real(1, :), 'b', t, Y_lin(1, :), 'g', t, Y_out, 'r', 'LineWidth', 1.5);
                
                legend({'model nieliniowy', ...
                        ['model liniowy' blanks(17) sprintf('E = %.3f', E_lin)], ...
                        ['model Hammersteina' blanks(5) sprintf('E = %.3f', E_out)]}, ...
                        'Location', 'best');
    
                title('Porównanie wyjścia obiektu testowego i jego modeli');
                ylabel('h [cm]');
                xlabel('t [s]');
                grid on;

            % --------------- LINEAR PH --------------- %
            else
                Y_stat = zeros(1, obiekt.kk);
                Y_fuzzy_Wa4 = evalfis(obj.linear_fis.Wa4, [U(1,:)' U(3,:)']);
                Y_fuzzy_Wb4 = evalfis(obj.linear_fis.Wb4, [U(1,:)' U(3,:)']);

                for k = 3:obiekt.kk
                    % Y_fuzzy_Wa4(k) = (obiekt.W_a1*(obiekt.Q_10+U(1,k-1)) + obiekt.W_a2*obiekt.Q_20 + obiekt.W_a3*(obiekt.Q_30+U(3,k-1)))./(obiekt.Q_10+obiekt.Q_20+obiekt.Q_30+U(1,k-1)+U(3,k-1));
                    % Y_fuzzy_Wb4(k) = (obiekt.W_b1*(obiekt.Q_10+U(1,k-1)) + obiekt.W_b2*obiekt.Q_20 + obiekt.W_b3*(obiekt.Q_30+U(3,k-1)))./(obiekt.Q_10+obiekt.Q_20+obiekt.Q_30+U(1,k-1)+U(3,k-1));
                    % Y_stat(k) = obiekt.pH_calc(Y_fuzzy_Wa4(k), Y_fuzzy_Wb4(k)) - obiekt.pH_0;
                    
                    Y_stat(k) = obiekt.pH_calc(Y_fuzzy_Wa4(k) + obiekt.W_a40, Y_fuzzy_Wb4(k) + obiekt.W_b40) - obiekt.pH_0;
                    Y_out(k) = -(a_pH.Q1 + a_pH.Q3)/2 * Y_out(k-1:-1:k-2)' + (-b_pH.Q1 + b_pH.Q3)/2 * Y_stat(k:-1:k-1)';
                end

                E_lin = sum((Y_real(2, :) - Y_lin(2, :)).^2) / obiekt.kk;
                E_out = sum((Y_real(2, :) - Y_out).^2) / obiekt.kk;
                fprintf("\nHAMMERSTEIN LINEAR MODEL\n");
                fprintf("E_lin = %.3f\n", E_lin);
                fprintf("E_out = %.3f\n", E_out);
    
                figure;
                plot(t, Y_real(2, :), 'b', t, Y_lin(2, :), 'g', t, Y_out, 'r', 'LineWidth', 1.5);

                legend({'model nieliniowy', ...
                        ['model liniowy' blanks(17) sprintf('E = %.3f', E_lin)], ...
                        ['model Hammersteina' blanks(5) sprintf('E = %.3f', E_out)]}, ...
                        'Location', 'best');

                title('Porównanie wyjścia obiektu testowego i jego modeli');
                ylabel('pH');
                xlabel('t [s]');
                grid on;
            end

            % saveas(gcf, sprintf('D:/EiTI/MGR/raporty/raport_MGR/pictures/HammersteinLinearModel_%d.png', index));  % Zapisuje jako plik PNG
        end

        function testNonlinearModel(obj, U, a_h, a_pH, b_h, b_pH, obiekt, type, index)
            
            [Y_real, Y_lin] = obiekt.modifiedEuler(U, obiekt.kk);
            t = 0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp; 
            Y_out = zeros(1, obiekt.kk);

            % --------------- NONLINEAR H --------------- %
            if(strcmp(type, 'h'))
                Y_fuzzy = evalfis(obj.nonlinear_fis.h, [sinh(U(1,:)'/15), sinh(U(3,:)'/15)]);
                for k = 2:obiekt.kk
                    Y_out(k) = - a_h*Y_out(k-1) + b_h*Y_fuzzy(k);
                end

                E_lin = sum((Y_real(1, :) - Y_lin(1, :)).^2) / obiekt.kk;
                E_out = sum((Y_real(1, :) - Y_out).^2) / obiekt.kk;
                fprintf("\nHAMMERSTEIN NONLINEAR MODEL\n");
                fprintf("E_lin = %.3f\n", E_lin);
                fprintf("E_out = %.3f\n", E_out);
    
                figure;
                plot(t, Y_real(1, :), 'b', t, Y_lin(1, :), 'g', t, Y_out, 'r', 'LineWidth', 1.5);
                
                legend({'model nieliniowy', ...
                        ['model liniowy' blanks(17) sprintf('E = %.3f', E_lin)], ...
                        ['model Hammersteina' blanks(5) sprintf('E = %.3f', E_out)]}, ...
                        'Location', 'best');
    
                title('Porównanie wyjścia obiektu testowego i jego modeli');
                ylabel('h [cm]');
                xlabel('t [s]');
                grid on;

            % --------------- NONLINEAR PH --------------- %
            else
                Y_stat = zeros(1, obiekt.kk);
                Y_fuzzy_Wa4 = evalfis(obj.nonlinear_fis.Wa4, [tanh(U(1,:)'/15) tanh(U(3,:)'/15)]);
                Y_fuzzy_Wb4 = evalfis(obj.nonlinear_fis.Wb4, [sinh(U(1,:)'/15) sinh(U(3,:)'/15)]);

                for k = 3:obiekt.kk
                    Y_stat(k) = obiekt.pH_calc(Y_fuzzy_Wa4(k) + obiekt.W_a40, Y_fuzzy_Wb4(k) + obiekt.W_b40) - obiekt.pH_0;
                    Y_out(k) = -(a_pH.Q1 + a_pH.Q3)/2 * Y_out(k-1:-1:k-2)' + (-b_pH.Q1 + b_pH.Q3)/2 * Y_stat(k:-1:k-1)';
                end

                E_lin = sum((Y_real(2, :) - Y_lin(2, :)).^2) / obiekt.kk;
                E_out = sum((Y_real(2, :) - Y_out).^2) / obiekt.kk;
                fprintf("\nHAMMERSTEIN NONLINEAR MODEL\n");
                fprintf("E_lin = %.3f\n", E_lin);
                fprintf("E_out = %.3f\n", E_out);
    
                figure;
                plot(t, Y_real(2, :), 'b', t, Y_lin(2, :), 'g', t, Y_out, 'r', 'LineWidth', 1.5);

                legend({'model nieliniowy', ...
                        ['model liniowy' blanks(17) sprintf('E = %.3f', E_lin)], ...
                        ['model Hammersteina' blanks(5) sprintf('E = %.3f', E_out)]}, ...
                        'Location', 'best');

                title('Porównanie wyjścia obiektu testowego i jego modeli');
                ylabel('pH');
                xlabel('t [s]');
                grid on;
                % saveas(gcf, sprintf('D:/EiTI/MGR/raporty/raport_MGR/pictures/HammersteinNonlinearModel_%d.png', index));  % Zapisuje jako plik PNG
            end
        end

        function checkStatic(obj, type)
            U = [linspace(obj.U_min, obj.U_max, 100);
                linspace(obj.U_min, obj.U_max, 100)];
            
            Y_Hout = zeros(100,100);
            Y_Aout = zeros(100,100);
            Y_Bout = zeros(100,100);

            if (strcmp(type, 'linear'))
                for i = 1:100
                    for j = 1:100
                        Y_Hout(i,j) = evalfis(obj.linear_fis.h, [U(1,j) U(2,i)]);
                        Y_Aout(i,j) = evalfis(obj.linear_fis.Wa4, [U(1,j) U(2,i)]);
                        Y_Bout(i,j) = evalfis(obj.linear_fis.Wb4, [U(1,j) U(2,i)]);
                    end
                end
            else
                for i = 1:100
                    for j = 1:100
                        Y_Hout(i,j) = evalfis(obj.nonlinear_fis.h, [sinh(U(1,j)/15) sinh(U(2,i)/15)]);
                        Y_Aout(i,j) = evalfis(obj.nonlinear_fis.Wa4, [tanh(U(1,j)/15) tanh(U(2,i)/15)]);
                        Y_Bout(i,j) = evalfis(obj.nonlinear_fis.Wb4, [sinh(U(1,j)/15) sinh(U(2,i)/15)]);
                    end
                end
            end

            figure;
            surf(obj.Q1_grid, obj.Q3_grid, Y_Hout);
            xlabel('Q_1 [ml/s]');
            ylabel('Q_3 [ml/s]');
            zlabel('h');
            title(sprintf('h'));
            shading interp;
            colorbar;
            
            figure;
            surf(obj.Q1_grid, obj.Q3_grid, Y_Aout);
            xlabel('Q_1 [ml/s]');
            ylabel('Q_3 [ml/s]');
            zlabel('pH');
            title(sprintf('Wa4'));
            shading interp;
            colorbar;
            
            figure;
            surf(obj.Q1_grid, obj.Q3_grid, Y_Bout);
            xlabel('Q_1 [ml/s]');
            ylabel('Q_3 [ml/s]');
            zlabel('pH');
            title(sprintf('Wb4'));
            shading interp;
            colorbar;
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