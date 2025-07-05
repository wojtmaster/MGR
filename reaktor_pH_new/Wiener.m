classdef Wiener < handle
    properties
        Y_h
        Y_Wa4
        Y_Wb4
        h
        Wa4
        Wb4
        pH
        Y_min = -30;
        Y_max = 30;
        U_min = -15;
        U_max= 15;

        linear_fis
        nonlinear_fis
    end

    methods

        function obj = Wiener(a_h, a_Wa4, a_Wb4, b_h, b_Wa4, b_Wb4, obiekt)
            U = [linspace(1.6, 31.6, 100)
                linspace(0.6, 30.6, 100)];
            [Q1_grid, Q3_grid] = meshgrid(U(1,:), U(2,:));
            obj.h = ((Q1_grid + obiekt.Q_20 + Q3_grid) / obiekt.C_V).^2 - obiekt.h_0;
            obj.Wa4 = (obiekt.W_a1*Q1_grid + obiekt.W_a2*obiekt.Q_20 + obiekt.W_a3*Q3_grid)./(Q1_grid+obiekt.Q_20+Q3_grid);
            obj.Wb4 = (obiekt.W_b1*Q1_grid + obiekt.W_b2*obiekt.Q_20 + obiekt.W_b3*Q3_grid)./(Q1_grid+obiekt.Q_20+Q3_grid);
            
            for i = 1:100
                for j = 1:100
                    obj.pH(i,j) = obiekt.pH_calc(obj.Wa4(i,j), obj.Wb4(i,j)) - obiekt.pH_0;
                end
            end

            obj.Wa4 = obj.Wa4 - obiekt.W_a40;
            obj.Wb4 = obj.Wb4 - obiekt.W_b40;

            U = [linspace(obj.U_min, obj.U_max, 100);
                linspace(obj.U_min, obj.U_max, 100)];

            [Q1_grid, Q3_grid] = meshgrid(U(1,:), U(2,:));

            obj.Y_h.Q1 = zeros(100, 100);
            obj.Y_h.Q3 = zeros(100, 100);
            obj.Y_Wa4.Q1 = zeros(100, 100);
            obj.Y_Wa4.Q3 = zeros(100, 100);
            obj.Y_Wb4.Q1 = zeros(100, 100);
            obj.Y_Wb4.Q3 = zeros(100, 100);
            
            y_h.Q1 = zeros(1, 100);
            y_h.Q3 = zeros(1, 100);
            y_Wa4.Q1 = zeros(1, 100);
            y_Wa4.Q3 = zeros(1, 100);
            y_Wb4.Q1 = zeros(1, 100);
            y_Wb4.Q3 = zeros(1, 100);
            
            for i = 1:100
                for j = 1:100
                    u = [ones(1,100) * U(1,j);
                        ones(1,100) * U(2,i)];
                    for k = 2:100
                        y_h.Q1(k) = - a_h*y_h.Q1(k-1) + b_h*u(1, k-1);
                        y_h.Q3(k) = - a_h*y_h.Q3(k-1) + b_h*u(2, k-1);

                        y_Wa4.Q1(k) = - a_Wa4*y_Wa4.Q1(k-1)' + b_Wa4*u(1, k-1)';
                        y_Wa4.Q3(k) = - a_Wa4*y_Wa4.Q3(k-1)' - b_Wa4*u(2, k-1)';
                    
                        y_Wb4.Q1(k) = - a_Wb4*y_Wb4.Q1(k-1)' + b_Wb4*u(1, k-1)';
                        y_Wb4.Q3(k) = - a_Wb4*y_Wb4.Q3(k-1)' + b_Wb4*u(2, k-1)';
                    end
                    obj.Y_h.Q1(i,j) = y_h.Q1(end);
                    obj.Y_h.Q3(i,j) = y_h.Q3(end);
                    obj.Y_Wa4.Q1(i,j) = y_Wa4.Q1(end);
                    obj.Y_Wa4.Q3(i,j) = y_Wa4.Q3(end);
                    obj.Y_Wb4.Q1(i,j) = y_Wb4.Q1(end);
                    obj.Y_Wb4.Q3(i,j) = y_Wb4.Q3(end);
                end
            end

            figure;
            surf(Q1_grid, Q3_grid, obj.Y_h.Q1, 'FaceColor', 'r'); % czerwona
            figure;
            surf(Q1_grid, Q3_grid, obj.Y_h.Q3, 'FaceColor', 'b');   % niebieska

            figure;
            surf(Q1_grid, Q3_grid, obj.Y_Wa4.Q1, 'FaceColor', 'r');
            figure;
            surf(Q1_grid, Q3_grid, obj.Y_Wa4.Q3, 'FaceColor', 'b');

            figure;
            surf(Q1_grid, Q3_grid, obj.Y_Wb4.Q1, 'FaceColor', 'r');
            figure;
            surf(Q1_grid, Q3_grid, obj.Y_Wb4.Q3, 'FaceColor', 'b');

        end

        function linearFuzzy(obj)
            % -------------------- h --------------------%
            X1 = obj.Y_h.Q1';  X2 = obj.Y_h.Q3';
            X = [X1(:) X2(:)];
            Y = obj.h';
            Y = Y(:);
            % Generowanie systemu rozmytego
            fis = genfis1([X Y], 5, 'gaussmf', 'linear');
            % Trening ANFIS
            options = anfisOptions('InitialFIS', fis, 'EpochNumber', 100, ...
                'DisplayANFISInformation', 0, 'DisplayErrorValues', 0);
            obj.linear_fis.h = anfis([X Y], options);

            % -------------------- pH --------------------%
            X1 = obj.Y_Wa4.Q1';  X2 = obj.Y_Wa4.Q3';
            X = [X1(:) X2(:)];
            Y = obj.Wa4';
            Y = Y(:);

            % Generowanie systemu rozmytego
            fis = genfis1([X Y], 5, 'gaussmf', 'linear');
            % Trening ANFIS
            options = anfisOptions('InitialFIS', fis, 'EpochNumber', 100, ...
                'DisplayANFISInformation', 0, 'DisplayErrorValues', 0);
            obj.linear_fis.Wa4 = anfis([X Y], options);

            X1 = obj.Y_Wb4.Q1';  X2 = obj.Y_Wb4.Q3';
            X = [X1(:) X2(:)];
            Y = obj.Wb4';
            Y = Y(:);

            % Generowanie systemu rozmytego
            fis = genfis1([X Y], 5, 'gaussmf', 'linear');
            % Trening ANFIS
            options = anfisOptions('InitialFIS', fis, 'EpochNumber', 100, ...
                'DisplayANFISInformation', 0, 'DisplayErrorValues', 0);
            obj.linear_fis.Wb4 = anfis([X Y], options);
        end

        function nonlinearFuzzy(obj)
            % ------------------ h ------------------ %
            X1 = obj.Y_h.Q1';  X2 = obj.Y_h.Q3';
            X = [sinh(X1(:)/15) sinh(X2(:)/15)];
            Y = obj.h';
            Y = Y(:);
            % Generowanie systemu rozmytego
            fis = genfis1([X Y], 3, 'gaussmf', 'linear');
            % Trening ANFIS
            options = anfisOptions('InitialFIS', fis, 'EpochNumber', 100, ...
                'DisplayANFISInformation', 0, 'DisplayErrorValues', 0);
            obj.nonlinear_fis.h = anfis([X Y], options);

            X1 = obj.Y_Wa4.Q1';  X2 = obj.Y_Wa4.Q3';
            X = [tanh(X1(:)/15) tanh(X2(:)/15)];
            Y = obj.Wa4';
            Y = Y(:);

            % Generowanie systemu rozmytego
            fis = genfis1([X Y], 3, 'gaussmf', 'linear');
            % Trening ANFIS
            options = anfisOptions('InitialFIS', fis, 'EpochNumber', 100, ...
                'DisplayANFISInformation', 0, 'DisplayErrorValues', 0);
            obj.nonlinear_fis.Wa4 = anfis([X Y], options);

            X1 = obj.Y_Wb4.Q1';  X2 = obj.Y_Wb4.Q3';
            X = [sinh(X1(:)/15) sinh(X2(:)/15)];
            Y = obj.Wb4';
            Y = Y(:);

            % Generowanie systemu rozmytego
            fis = genfis1([X Y], 3, 'gaussmf', 'linear');
            % Trening ANFIS
            options = anfisOptions('InitialFIS', fis, 'EpochNumber', 100, ...
                'DisplayANFISInformation', 0, 'DisplayErrorValues', 0);
            obj.nonlinear_fis.Wb4 = anfis([X Y], options);
        end

        function testLinearModel(obj, U, a_h, a_Wa4, a_Wb4, b_h, b_Wa4, b_Wb4, obiekt, type, index)
            [Y_real, Y_lin] = obiekt.modifiedEuler(U, obiekt.kk); % Symulacja rzeczywistego układu
            t = 0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp; 

            Y_out.Q1 = zeros(1, obiekt.kk);
            Y_out.Q3 = zeros(1, obiekt.kk);
            Y_fuzzy = zeros(1, obiekt.kk);

            % -------------------- LINEAR H -------------------- %
            if(strcmp(type, 'h'))
                for k = 2:obiekt.kk
                    Y_out.Q1(k) = - a_h*Y_out.Q1(k-1) + b_h * U(1,k-1);
                    Y_out.Q3(k) = - a_h*Y_out.Q3(k-1) + b_h * U(3,k-1);

                    Y_fuzzy(k) = ((obiekt.Q_10+obiekt.Q_20+obiekt.Q_30+Y_out.Q1(k)+Y_out.Q3(k)) / obiekt.C_V).^2 - obiekt.h_0;
                    % Y_fuzzy(k) = evalfis(obj.linear_fis.h, [Y_out.Q1(k), Y_out.Q3(k)]);
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
                y_Wa4.Q1 = zeros(1,obiekt.kk);
                y_Wa4.Q3 = zeros(1,obiekt.kk);
                y_Wb4.Q1 = zeros(1,obiekt.kk);
                y_Wb4.Q3 = zeros(1,obiekt.kk);
                for k = 4:obiekt.kk
                    y_Wa4.Q1(k) = -a_Wa4 * y_Wa4.Q1(k-1)' + b_Wa4*U(1,k-1)';
                    y_Wa4.Q3(k) = -a_Wa4 * y_Wa4.Q3(k-1)'- b_Wa4*U(3,k-1)';
                    y_Wb4.Q1(k) = -a_Wb4 * y_Wb4.Q1(k-1)' + b_Wb4*U(1,k-1)';
                    y_Wb4.Q3(k) = -a_Wb4 * y_Wb4.Q3(k-1)' + b_Wb4*U(3,k-1)';

                    % W_a4 = (obiekt.W_a1*(obiekt.Q_10+y_Wa4.Q1(k))+ obiekt.W_a2*obiekt.Q_20+obiekt.W_a3*(obiekt.Q_30-y_Wa4.Q3(k)))./(obiekt.Q_10+obiekt.Q_20+obiekt.Q_30+y_Wa4.Q1(k)-y_Wa4.Q3(k));
                    % W_b4 = (obiekt.W_b1*(obiekt.Q_10-y_Wb4.Q1(k)) + obiekt.W_b2*obiekt.Q_20 + obiekt.W_b3*(obiekt.Q_30-y_Wb4.Q3(k)))./(obiekt.Q_10+obiekt.Q_20+obiekt.Q_30-y_Wb4.Q1(k)-y_Wb4.Q3(k));
                    % Y_fuzzy(k) = obiekt.pH_calc(W_a4, W_b4) - obiekt.pH_0;

                    W_a4 = evalfis(obj.linear_fis.Wa4, [y_Wa4.Q1(k), y_Wa4.Q3(k)]);
                    W_b4 = evalfis(obj.linear_fis.Wb4, [y_Wb4.Q1(k), y_Wb4.Q3(k)]);
                    Y_fuzzy(k) = obiekt.pH_calc(W_a4 + obiekt.W_a40, W_b4 + obiekt.W_b40) - obiekt.pH_0;
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

        function testNonlinearModel(obj, U, a_h, a_Wa4, a_Wb4, b_h, b_Wa4, b_Wb4, obiekt, type, index)
            [Y_real, Y_lin] = obiekt.modifiedEuler(U, obiekt.kk); % Symulacja rzeczywistego układu
            t = 0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp;

            Y_out.Q1 = zeros(1,obiekt.kk);
            Y_out.Q3 = zeros(1,obiekt.kk);
            Y_fuzzy = zeros(1,obiekt.kk);

            % -------------------- NONLINEAR H -------------------- %
            if(strcmp(type, 'h'))
                for k = 2:obiekt.kk
                    Y_out.Q1(k) = - a_h*Y_out.Q1(k-1) + b_h * U(1,k-1);
                    Y_out.Q3(k) = - a_h*Y_out.Q3(k-1) + b_h * U(3,k-1);
                    Y_fuzzy(k) = evalfis(obj.nonlinear_fis.h, [sinh(Y_out.Q1(k)/15), sinh(Y_out.Q3(k)/15)]);
                end
                E_lin = sum((Y_real(1, :) - Y_lin(1, :)).^2) / obiekt.kk;
                E_out = sum((Y_real(1, :) - Y_fuzzy).^2) / obiekt.kk;
                fprintf("\nWIENER NONLINEAR MODEL\n");
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
                y_Wa4.Q1 = zeros(1,obiekt.kk);
                y_Wa4.Q3 = zeros(1,obiekt.kk);
                y_Wb4.Q1 = zeros(1,obiekt.kk);
                y_Wb4.Q3 = zeros(1,obiekt.kk);
                for k = 4:obiekt.kk
                    y_Wa4.Q1(k) = -a_Wa4 * y_Wa4.Q1(k-1)' + b_Wa4*U(1,k-1)';
                    y_Wa4.Q3(k) = -a_Wa4 * y_Wa4.Q3(k-1)'- b_Wa4*U(3,k-1)';
                    y_Wb4.Q1(k) = -a_Wb4 * y_Wb4.Q1(k-1)' + b_Wb4*U(1,k-1)';
                    y_Wb4.Q3(k) = -a_Wb4 * y_Wb4.Q3(k-1)' + b_Wb4*U(3,k-1)';

                    W_a4 = evalfis(obj.nonlinear_fis.Wa4, [tanh(y_Wa4.Q1(k)/15), tanh(y_Wa4.Q3(k)/15)]);
                    W_b4 = evalfis(obj.nonlinear_fis.Wb4, [sinh(y_Wb4.Q1(k)/15), sinh(y_Wb4.Q3(k)/15)]);
                
                    Y_fuzzy(k) = obiekt.pH_calc(W_a4 + obiekt.W_a40, W_b4 + obiekt.W_b40) - obiekt.pH_0;
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

        function checkStatic(obj, type)
            U = [linspace(obj.U_min, obj.U_max, 100);
                linspace(obj.U_min, obj.U_max, 100)];
            [Q1_grid, Q3_grid] = meshgrid(U(1,:), U(2,:));
            
            Y_Hout = zeros(100,100);
            Y_Aout = zeros(100,100);
            Y_Bout = zeros(100,100);

            if (strcmp(type, 'linear'))
                for i = 1:100
                    for j = 1:100
                        Y_Hout(i,j) = evalfis(obj.linear_fis.h, [obj.Y_h.Q1(i,j), obj.Y_h.Q3(i,j)]);
                        Y_Aout(i,j) = evalfis(obj.linear_fis.Wa4, [obj.Y_Wa4.Q1(i,j), obj.Y_Wa4.Q3(i,j)]);
                        Y_Bout(i,j) = evalfis(obj.linear_fis.Wb4, [obj.Y_Wb4.Q1(i,j), obj.Y_Wb4.Q3(i,j)]);
                    end
                end
            else
                for i = 1:100
                    for j = 1:100
                        Y_Hout(i,j) = evalfis(obj.nonlinear_fis.h, [sinh(obj.Y_h.Q1(i,j)/15), sinh(obj.Y_h.Q3(i,j)/15)]);
                        Y_Aout(i,j) = evalfis(obj.nonlinear_fis.Wa4, [tanh(obj.Y_Wa4.Q1(i,j)/15), tanh(obj.Y_Wa4.Q3(i,j)/15)]);
                        Y_Bout(i,j) = evalfis(obj.nonlinear_fis.Wb4, [sinh(obj.Y_Wb4.Q1(i,j)/15), sinh(obj.Y_Wb4.Q3(i,j)/15)]);
                    end
                end
            end

            figure;
            surf(Q1_grid, Q3_grid, Y_Hout);
            xlabel('Q_1 [ml/s]');
            ylabel('Q_3 [ml/s]');
            zlabel('h');
            title(sprintf('h'));
            shading interp;
            colorbar;
            
            figure;
            surf(Q1_grid, Q3_grid, Y_Aout);
            xlabel('Q_1 [ml/s]');
            ylabel('Q_3 [ml/s]');
            zlabel('pH');
            title(sprintf('Wa4'));
            shading interp;
            colorbar;
            
            figure;
            surf(Q1_grid, Q3_grid, Y_Bout);
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
            title(sprintf('Funkcje przynależności dla Y - następniki %s', s));
            xlabel('y');
            ylabel('$\mu(y)$', 'Interpreter', 'latex');
            grid on;

            % saveas(gcf, sprintf('D:/EiTI/MGR/raporty/raport_MGR/pictures/WienerfuzzySets_%s.png', s));  % Zapisuje jako plik PNG
        end

    end
end
