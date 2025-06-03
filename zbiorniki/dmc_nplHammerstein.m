function [y_mod, u, E, E_u, E_y] = dmc_nplHammerstein(obj, y_zad, u_D, fis, delay, kk, F_10, F_D0, obiekt, type)
            % Alokacja pamięci
            u = zeros(2, kk);
            u(2,:) = u_D;

            delta_up = zeros(1, obj.D-1);
            delta_uz = zeros(1, obj.D_disturbance);
            delta_uk = 0;
            delta_u = zeros(1, kk);
            
            % Sterowanie DMC
            y_mod = zeros(size(y_zad));
            y_mod(1:delay+1) = y_zad(1:delay+1);
            y_fuzzy = zeros(size(y_zad));

            h_20 = obiekt.h_20;
            V_1 = zeros(size(y_zad));
            V_2 = zeros(size(y_zad));
            V_1(1:delay+1) = obiekt.V_10;
            V_2(1:delay+1) = obiekt.V_20;
            V_1e = zeros(size(y_zad));
            V_2e = zeros(size(y_zad));
            V_1e(1:delay+1) = obiekt.V_10;
            V_2e(1:delay+1) = obiekt.V_20;

            % Główna pętla sterowania
            for k = delay+2:kk
                Symulacja nieliniowego modelu procesu
                [~, V_1(k), V_2(k), V_1e(k), V_2e(k)] = obiekt.realSimulation(...
                    [u(1, k-(delay+1)) u(2, k-1)]', F_10, F_D0, h_20, ...
                    V_1(k-1), V_2(k-1), V_1e(k-1), V_2e(k-1));
                
                Obliczenie wzmocnienia z modelu Hammersteina
                if (strcmp(type, 'linear'))
                    y_fuzzy(k) = evalfis(fis, [u(1,k-(delay+1)), u(2,k-1)]);
                else
                    gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));
        
                    deg_u1 = [gaussmf_val(u(1,k-(delay+1)), fis.sigma_F1, fis.F1_center(1)), ...
                             gaussmf_val(u(1,k-(delay+1)), fis.sigma_F1, fis.F1_center(2)), ...
                             gaussmf_val(u(1,k-(delay+1)), fis.sigma_F1, fis.F1_center(3))];
                    deg_u2 = [gaussmf_val(u(2,k-1), fis.sigma_FD, fis.FD_center(1)), ...
                             gaussmf_val(u(2,k-1), fis.sigma_FD, fis.FD_center(2))];
                    
                    degrees_all(1) = deg_u1(1) * deg_u2(1);
                    degrees_all(2) = deg_u1(1) * deg_u2(2);
                    degrees_all(3) = deg_u1(2) * deg_u2(1);
                    degrees_all(4) = deg_u1(2) * deg_u2(2);
                    degrees_all(5) = deg_u1(3) * deg_u2(1);
                    degrees_all(6) = deg_u1(3) * deg_u2(2);
                    
                    w = 0;
                    output = 0;
                    for i = 1:fis.rules_number
                        output = output + degrees_all(i)*(fis.a_param(i)*sinh(u(1,k-(delay+1))/22.5) + ...
                                fis.b_param(i)*sinh(u(2,k-1)/7.5) + fis.c_param(i));
                        w = w + degrees_all(i);
                    end
                    y_fuzzy(k) = output / w;
                end
        
                Obliczenie wzmocnienia statycznego
                if (abs(u(1,k-(delay+1)) + u(2,k-1)) >= 1e-3)
                    gain = y_fuzzy(k) / (u(1,k-(delay+1)) + u(2,k-1));
                else
                    gain = 1; % Wartość domyślna dla małych sygnałów
                end
                
                Linearyzacja w punkcie pracy i aktualizacja macierzy DMC
                obiekt.linearization(F_10 + u(1, k-(delay+1)), F_D0 + u(2, k-1));
                [~, ~, s, ~] = obiekt.mse();
                
                Aktualizacja macierzy dynamiki
                obj.dynamic_matrix(s');
                
                Predykcja swobodnej odpowiedzi
                Y_0 = obiekt.freeAnswer([u(1, k-(delay+1)) u(2, k-1)]', F_10, F_D0, h_20, ...
                                      V_1(k), V_2(k), V_1e(k), V_2e(k), obj.N);

                Y_zad = y_zad(k) * ones(obj.N, 1);
                delta_uk = obj.K * (Y_zad - Y_0);
                
                for i = 1:obj.Nu
                    Ograniczenia przyrostów sterowania
                    delta_uk = min(delta_uk, obj.delta_uk_max);
                    delta_uk = max(delta_uk, -obj.delta_uk_max);
                end

                Aktualizacja sterowania
                u(1, k) = u(1, k-1) + delta_uk(1);

                Ograniczenia sterowania
                u(1, k) = min(u(1, k), obj.u_max);
                u(1, k) = max(u(1, k), -obj.u_max);
                
                Zapamiętanie przyrostu sterowania
                delta_u(k) = delta_uk(1);
                
                Predykcja wyjścia
                Y = Y_0 + gain * obj.M * delta_uk;

                disp(u(1, k-(delay+1)));
                disp(u(2, k-1));
                if (abs(sum(Y)) > 0.1 && abs(sum(Y_0)) > 0.1)
                    figure;
                    plot(Y_0, 'b', 'LineWidth', 2);
                    hold on;
                    plot(obj.M * delta_uk, 'r')
                    plot(Y, 'g')
                    keyboard;
                    close all;
                end

                y_mod(k) = Y(1);
                
                Ograniczenia wyjścia
                y_mod(k) = min(y_mod(k), obj.y_max);
                y_mod(k) = max(y_mod(k), -obj.y_max);
            end
            E_u = obj.lambda .* sum(delta_u.^2)/kk;
            E_y = sum((y_zad - y_mod).^2)/kk;
            E =  E_u + E_y;
        end