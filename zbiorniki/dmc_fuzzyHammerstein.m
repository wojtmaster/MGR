        function [y_mod, u, E, E_u, E_y] = dmc_fuzzyHammerstein(obj, y_zad, u_D, a, b, fis, delay, kk, F_10, F_D0, h_20, obiekt, type)
            gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));
            %% Alokacja pamięci
            u = zeros(2, kk);
            u(2,:) = u_D;

            delta_up = zeros(1, obj.D-1);
            delta_uz = zeros(1, obj.D_disturbance);
            delta_uk = 0;
            delta_u = zeros(1, kk);

            F_1 = [-30 0 30];
            F_D = [-10 0 10];
            L = [1 1 1];
            sigma = 15;
            rules_number = 3;
            
            Y = zeros(1, rules_number);
            KE = zeros(1, rules_number);
            KU = cell(1, rules_number);
            KZ = cell(1, rules_number);

            for i = 1:rules_number
                obj.lambda = L(i);
                obiekt.linearization(F_10+F_1(i), F_D0+F_D(i));
                [~, ~, ~, ~, s, s_disturbance] = obiekt.mse();
                obj.dynamic_matrix(s');
                obj.past_matrix(s');
                obj.matrix_disturbance(s_disturbance');
        
                Y(i) = obiekt.h_20 - h_20;
                KE(i) = obj.ke;
                KU{i} = obj.ku;
                KZ{i} = obj.kz;
            end

            %% Sterowanie DMC
            y_fuzzy = zeros(size(y_zad));
            y_mod = zeros(size(y_zad));
            y_mod(1:delay+2) = y_zad(1:delay+2);
            for k = delay+2:kk
                if (strcmp(type, 'linear'))
                    y_fuzzy(k) = evalfis(fis, [u(1,k-(delay+1)), u(2,k-1)]);
                else
                    gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));

                    deg_u1 = [gaussmf_val(u(1,k-(delay+1)), fis.sigma_F1, fis.F1_center(1)), gaussmf_val(u(1,k-(delay+1)), fis.sigma_F1, fis.F1_center(2)), gaussmf_val(u(1,k-(delay+1)), fis.sigma_F1, fis.F1_center(3))];
                    deg_u2 = [gaussmf_val(u(2,k-1), fis.sigma_FD, fis.FD_center(1)), gaussmf_val(u(2,k-1), fis.sigma_FD, fis.FD_center(2))];
                    
                    degrees_all(1) = deg_u1(1) * deg_u2(1);
                    degrees_all(2) = deg_u1(1) * deg_u2(2);
                    degrees_all(3) = deg_u1(2) * deg_u2(1);
                    degrees_all(4) = deg_u1(2) * deg_u2(2);
                    degrees_all(5) = deg_u1(3) * deg_u2(1);
                    degrees_all(6) = deg_u1(3) * deg_u2(2);
                    
                    w = 0;
                    output = 0;
                    for i = 1:fis.rules_number
                        output = output + degrees_all(i)*(fis.a_param(i)*sinh(u(1,k-(delay+1))/22.5) + fis.b_param(i)*sinh(u(2,k-1)/7.5) + fis.c_param(i));
                        w = w + degrees_all(i);
                    end
                    y_fuzzy(k) = output / w;
                end
                
                if(abs(u(1,k-(delay+1)) + u(2,k-1)) > 1e-3)
                    gain = y_fuzzy(k) / (u(1,k-(delay+1)) + u(2,k-1));
                else
                    gain = 1;
                end

                y_mod(k) = - a*y_mod(k-1:-1:k-2)' + gain*b*u(1, k-(delay+1))' + gain*b*u(2, k-1)';
                % Ograniczenia wartości sygnału wyjściowego, tj. wysokości h_2
                y_mod(k) = min(y_mod(k), obj.y_max);
                y_mod(k) = max(y_mod(k), -obj.y_max);

                % Przepisanie sterowań do wektora przeszłych sterowań
                delta_up = [delta_uk, delta_up(1:end-1)];
                delta_uz = [u(2, k) - u(2, k-1) , delta_uz(1:end-1)];

                % Oblicznie uchybu    
                e = y_zad(k) - y_mod(k);

                output = 0;
                w = 0;
                % strength = strength / sum(strength);
                for i = 1:rules_number
                    degree = gaussmf_val(y_mod(k), sigma, Y(i));
                    output = output + degree * (KE(i) * e - KU{i} * delta_up' - KZ{i} * delta_uz');
                    w = w + degree;
                end
                delta_uk = output / w;

                % Ograniczenie wartości przyrostu sterowania
                delta_uk = min(delta_uk, obj.delta_uk_max);
                delta_uk = max(delta_uk, -obj.delta_uk_max);
                
                % Obliczenie sterowania dla chwili (i+1) w chwili i-tej
                u(1, k) = u(1, k-1) + delta_uk;
                
                % Ograniczenie sterowania
                u(1, k) = min(u(1, k), obj.u_max);
                u(1, k) = max(u(1, k), -obj.u_max);

                delta_u(k) = delta_uk;
            end
            E_u = obj.lambda .* sum(delta_u.^2)/kk;
            E_y = sum((y_zad - y_mod).^2)/kk;
            E =  E_u + E_y;
        end