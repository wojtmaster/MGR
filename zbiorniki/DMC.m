classdef DMC < handle
    properties
        N
        Nu
        D
        D_disturbance
        M
        M_p
        M_zP
        K
        ke
        ku
        kz
        lambda

        u_max = 45;
        delta_uk_max = 2.5;
        y_max = 25;
    end

    methods
        function obj = DMC(N, Nu, D, D_disturbance, lambda)
            obj.N = N;
            obj.Nu = Nu;
            obj.D = D;
            obj.D_disturbance = D_disturbance;
            obj.lambda = lambda;
        end

        function [] = dynamic_matrix(obj, s)
            obj.M = zeros(obj.N, obj.Nu);
            
            % Implementacja macierzy M
            for i = 1:obj.N
                for j = 1:obj.Nu
                    if(i >= j)
                        obj.M(i,j) = s(i-j+1);
                    end
                end
            end

            % Macierz K
            obj.K = ((obj.M' * obj.M + obj.lambda * eye(obj.Nu))^(-1)) * obj.M';
            obj.ke = sum(obj.K(1, :));          
        end

        function [] = past_matrix(obj, s)
            obj.M_p = zeros(obj.N, obj.D-1);
        
            % Implementacja macierzy M_p
            for i = 1:obj.N
                for j = 1:obj.D-1
                    if(i+j <= obj.D)
                        obj.M_p(i,j) = s(i+j) - s(j);
                    else 
                        obj.M_p(i,j) = s(obj.D) - s(j);
                    end
                end
            end
            obj.ku = obj.K(1, :) * obj.M_p;
        end

        function [] = matrix_disturbance(obj, s)
            obj.M_zP = zeros(obj.N, obj.D_disturbance-1);

            % Implementacja macierzy M_zP
            for i = 1:obj.N
                for j = 1:obj.D_disturbance-1
                    if(i+j <= obj.D_disturbance)
                        obj.M_zP(i,j) = s(i+j) - s(j);
                    else 
                        obj.M_zP(i,j) = s(obj.D_disturbance) - s(j);
                    end
                end
            end

            obj.M_zP = [s(1:obj.N), obj.M_zP];
            obj.kz = obj.K(1, :) * obj.M_zP;
        end
    
        function [y_mod, u, E, E_u, E_y] = dmc_analiticHammerstein(obj, y_zad, u_D, a, b, fis, delay, kk, type)
            %% Alokacja pamięci
            u = zeros(2, kk);
            u(2,:) = u_D;

            delta_up = zeros(1, obj.D-1);
            delta_uz = zeros(1, obj.D_disturbance);
            delta_uk = 0;
            delta_u = zeros(1, kk);
            
            %% Sterowanie DMC
            y_fuzzy = zeros(size(y_zad));
            y_mod = zeros(size(y_zad));
            y_mod(1:delay+1) = y_zad(1:delay+1);

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
                
                if(u(1,k-(delay+1)) ~= 0)
                    gain = y_fuzzy(k) / (u(1,k-(delay+1)) + u(2,k-1));
                else
                    gain = 1;
                end

                y_mod(k) = - a*y_mod(k-1:-1:k-2)' + gain * b *u(1, k-(delay+1)) + gain * b*u(2, k-1);

                % Ograniczenia wartości sygnału wyjściowego, tj. wysokości h_2
                y_mod(k) = min(y_mod(k), obj.y_max);
                y_mod(k) = max(y_mod(k), -obj.y_max);

                % Przepisanie sterowań do wektora przeszłych sterowań
                delta_up = [delta_uk, delta_up(1:end-1)];
                delta_uz = [u(2, k) - u(2, k-1) , delta_uz(1:end-1)];

                % Oblicznie uchybu    
                e = y_zad(k) - y_mod(k);
            
                % Obliczenie przyrostu sterowania dla chwili (i+1) w chwili i-tej
                delta_uk = obj.ke * e - obj.ku * delta_up' - obj.kz * delta_uz';
                
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

        function [y_mod, u, E, E_u, E_y] = dmc_numericHammerstein(obj, y_zad, u_D, a, b, fis, delay, kk, type)
            % Alokacja pamięci
            u = zeros(2, kk);
            u(2,:) = u_D;
            
            delta_up = zeros(obj.D-1,1);
            delta_uk = zeros(obj.Nu,1);
            delta_uz = zeros(obj.D_disturbance, 1);
            delta_u = zeros(1, kk);

            Y_max = ones(obj.N,1)*obj.y_max;
            Y_min = ones(obj.N,1)*(-obj.y_max);
            
            U_max = ones(obj.Nu,1)*obj.u_max;
            U_min = ones(obj.Nu,1)*(-obj.u_max);
            delta_U_max = ones(obj.Nu,1)*obj.delta_uk_max;
            delta_U_min = ones(obj.Nu,1)*(-obj.delta_uk_max);
            
            J = tril(ones(obj.Nu));
            A = [-J; J; -obj.M; obj.M];
            H = 2*(obj.M'*obj.M + obj.lambda*eye(obj.Nu));
            
            y_fuzzy = zeros(size(y_zad));
            y_mod = zeros(size(y_zad));
            y_mod(1:delay+1) = y_zad(1:delay+1);

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

                if(u(1,k-(delay+1)) ~= 0)
                    gain = y_fuzzy(k) / (u(1,k-(delay+1)) + u(2,k-1));
                else
                    gain = 1;
                end

                y_mod(k) = - a*y_mod(k-1:-1:k-2)' + gain*b*u(1, k-(delay+1)) + gain*b*u(2, k-1);
                   
                % Przepisanie sterowań do wektora przeszłych sterowań
                delta_up = [delta_uk(1); delta_up(1:end-1)];
                delta_uz = [u(2, k) - u(2, k-1); delta_uz(1:end-1)];
            
                U_k_1 = ones(obj.Nu,1)*u(1, k-1);
                Y = ones(obj.N,1)*y_mod(k);
                Y_0 = Y + obj.M_p*delta_up + obj.M_zP*delta_uz;
                Y_zad = ones(obj.N,1)*y_zad(k);
                B = [-U_min + U_k_1;
                    U_max - U_k_1;
                    -Y_min + Y_0;
                    Y_max - Y_0];
                f = -2*obj.M'*(Y_zad - Y_0);
                
                % Obliczenie przyrostu sterowania dla chwili (i+1) w chwili i-tej
                options = optimoptions('quadprog', 'Display', 'off');
                delta_uk = quadprog(H,f,A,B,[],[],delta_U_min,delta_U_max,[],options);
                
                % Obliczenie sterowania
                u(1,k) = u(1,k-1) + delta_uk(1);
                delta_u(k) = delta_uk(1);
            end
            E_u = obj.lambda .* sum(delta_u.^2)/kk;
            E_y = sum((y_zad - y_mod).^2)/kk;
            E =  E_u + E_y;
        end

        function [y_mod, u, E, E_u, E_y] = dmc_slHammerstein(obj, y_zad, u_D, fis, delay, kk, F_10, F_D0, obiekt, type)
            %% Alokacja pamięci
            u = zeros(2, kk);
            u(2,:) = u_D;

            delta_up = zeros(1, obj.D-1);
            delta_uz = zeros(1, obj.D_disturbance);
            delta_uk = 0;
            delta_u = zeros(1, kk);
            
            %% Sterowanie DMC
            y_mod = zeros(size(y_zad));
            y_mod(1:delay+1) = y_zad(1:delay+1);
            y_fuzzy = zeros(size(y_zad));
            for k = delay+2:kk
                obiekt.linearization(F_10 + u(1,k-(delay+1)), F_D0 + u(2,k-1));
                [a, b, s, s_disturbance] = obiekt.mse();
                obj.dynamic_matrix(s');
                obj.past_matrix(s');
                obj.matrix_disturbance(s_disturbance');

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

                if(u(1,k-(delay+1)) + u(2,k-1) ~= 0)
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
            
                % Obliczenie przyrostu sterowania dla chwili (i+1) w chwili i-tej
                delta_uk = obj.ke * e - obj.ku * delta_up' - obj.kz * delta_uz';
                
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
                % Symulacja nieliniowego modelu procesu
                [~, V_1(k), V_2(k), V_1e(k), V_2e(k)] = obiekt.realSimulation(...
                    [u(1, k-(delay+1)) u(2, k-1)]', F_10, F_D0, h_20, ...
                    V_1(k-1), V_2(k-1), V_1e(k-1), V_2e(k-1));
                
                % Obliczenie wzmocnienia z modelu Hammersteina
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
        
                % Obliczenie wzmocnienia statycznego
                if (abs(u(1,k-(delay+1)) + u(2,k-1)) >= 1e-3)
                    gain = y_fuzzy(k) / (u(1,k-(delay+1)) + u(2,k-1));
                else
                    gain = 1; % Wartość domyślna dla małych sygnałów
                end
                
                % Linearyzacja w punkcie pracy i aktualizacja macierzy DMC
                obiekt.linearization(F_10 + u(1, k-(delay+1)), F_D0 + u(2, k-1));
                [~, ~, s, ~] = obiekt.mse();
                
                % Aktualizacja macierzy dynamiki
                obj.dynamic_matrix(s');
                
                % Predykcja swobodnej odpowiedzi
                Y_0 = obiekt.freeAnswer([u(1, k-(delay+1)) u(2, k-1)]', F_10, F_D0, h_20, ...
                                      V_1(k), V_2(k), V_1e(k), V_2e(k), obj.N);

                Y_zad = y_zad(k) * ones(obj.N, 1);
                delta_uk = obj.K * (Y_zad - Y_0);
                
                for i = 1:obj.Nu
                    % Ograniczenia przyrostów sterowania
                    delta_uk = min(delta_uk, obj.delta_uk_max);
                    delta_uk = max(delta_uk, -obj.delta_uk_max);
                end

                % Aktualizacja sterowania
                u(1, k) = u(1, k-1) + delta_uk(1);

                % Ograniczenia sterowania
                u(1, k) = min(u(1, k), obj.u_max);
                u(1, k) = max(u(1, k), -obj.u_max);
                
                % Zapamiętanie przyrostu sterowania
                delta_u(k) = delta_uk(1);
                
                % Predykcja wyjścia
                Y = Y_0 + gain * obj.M * delta_uk;

                % disp(u(1, k-(delay+1)));
                % disp(u(2, k-1));
                % if (abs(sum(Y)) > 0.1 && abs(sum(Y_0)) > 0.1)
                %     figure;
                %     plot(Y_0, 'b', 'LineWidth', 2);
                %     hold on;
                %     plot(obj.M * delta_uk, 'r')
                %     plot(Y, 'g')
                %     keyboard;
                %     close all;
                % end

                y_mod(k) = Y(1);
                
                % Ograniczenia wyjścia
                y_mod(k) = min(y_mod(k), obj.y_max);
                y_mod(k) = max(y_mod(k), -obj.y_max);
            end
            E_u = obj.lambda .* sum(delta_u.^2)/kk;
            E_y = sum((y_zad - y_mod).^2)/kk;
            E =  E_u + E_y;
        end

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
            L = [10 0.5 0.5];
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

        function [y_out, u, E, E_u, E_y] = dmc_analiticWiener(obj, y_zad, u_D, a, b, fis, delay, kk, type)
            %% Alokacja pamięci
            u = zeros(2, kk);
            u(2,:) = u_D;

            delta_up = zeros(1, obj.D-1);
            delta_uz = zeros(1, obj.D_disturbance);
            delta_uk = 0;
            delta_u = zeros(1, kk);
            
            %% Sterowanie DMC
            y_mod = zeros(size(y_zad));
            y_out = zeros(size(y_zad));
            y_mod(1:delay+1) = y_zad(1:delay+1);
            for k = delay+2:kk
                y_mod(k) = - a*y_mod(k-1:-1:k-2)' + b*u(1, k-(delay+1)) + b*u(2, k-1);

                if (strcmp(type, 'linear'))
                    y_out(k) = evalfis(fis, y_mod(k));
                else
                    gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));

                    output = 0;
                    w = 0;
                    degrees = [gaussmf_val(y_mod(k), fis.sigma, fis.Y_center(1)), ...
                               gaussmf_val(y_mod(k), fis.sigma, fis.Y_center(2)), ...
                               gaussmf_val(y_mod(k), fis.sigma, fis.Y_center(3))];
                    for i = 1:fis.rules_number
                        output = output +  degrees(i) * (fis.a_param(i)*sinh(y_mod(k)/36) + fis.b_param(i));
                        w = w + degrees(i);
                    end
                    y_out(k) = output / w;
                end

                % Ograniczenia wartości sygnału wyjściowego, tj. wysokości h_2
                y_out(k) = min(y_out(k), obj.y_max);
                y_out(k) = max(y_out(k), -obj.y_max);

                % Przepisanie sterowań do wektora przeszłych sterowań
                delta_up = [delta_uk, delta_up(1:end-1)];
                delta_uz = [u(2, k) - u(2, k-1) , delta_uz(1:end-1)];

                % Oblicznie uchybu    
                e = y_zad(k) - y_out(k);
            
                % Obliczenie przyrostu sterowania dla chwili (i+1) w chwili i-tej
                delta_uk = obj.ke * e - obj.ku * delta_up' - obj.kz * delta_uz';
                
                % Ograniczenie wartości przyrostu sterowania
                delta_uk = min(delta_uk, obj.delta_uk_max);
                delta_uk = max(delta_uk, -obj.delta_uk_max);
                delta_u(k) = delta_uk;
                
                % Obliczenie sterowania dla chwili (i+1) w chwili i-tej
                u(1, k) = u(1, k-1) + delta_uk;
                
                % Ograniczenie sterowania
                u(1, k) = min(u(1, k), obj.u_max);
                u(1, k) = max(u(1, k), -obj.u_max);
            end
            E_u = obj.lambda .* sum(delta_u.^2)/kk;
            E_y = sum((y_zad - y_out).^2)/kk;
            E =  E_u + E_y;
        end

        function [y_out, u, E, E_u, E_y] = dmc_numericWiener(obj, y_zad, u_D, a, b, fis, delay, kk, type)
            % Alokacja pamięci
            u = zeros(2, kk);
            u(2,:) = u_D;
            
            delta_up = zeros(obj.D-1,1);
            delta_uk = zeros(obj.Nu,1);
            delta_uz = zeros(obj.D_disturbance, 1);
            delta_u = zeros(1, kk);

            Y_max = ones(obj.N,1)*obj.y_max;
            Y_min = ones(obj.N,1)*(-obj.y_max);
            
            U_max = ones(obj.Nu,1)*obj.u_max;
            U_min = ones(obj.Nu,1)*(-obj.u_max);
            delta_U_max = ones(obj.Nu,1)*obj.delta_uk_max;
            delta_U_min = ones(obj.Nu,1)*(-obj.delta_uk_max);
            
            J = tril(ones(obj.Nu));
            A = [-J; J; -obj.M; obj.M];
            H = 2*(obj.M'*obj.M + obj.lambda*eye(obj.Nu));
            
            y_mod = zeros(size(y_zad));
            y_out = zeros(size(y_zad));
            y_mod(1:delay+1) = y_zad(1:delay+1);
            for k = delay+2:kk
                y_mod(k) = - a*y_mod(k-1:-1:k-2)' + b*u(1, k-(delay+1)) + b*u(2, k-1);

                if (strcmp(type, 'linear'))
                    y_out(k) = evalfis(fis, y_mod(k));
                else
                    gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));

                    output = 0;
                    w = 0;
                    degrees = [gaussmf_val(y_mod(k), fis.sigma, fis.Y_center(1)), ...
                               gaussmf_val(y_mod(k), fis.sigma, fis.Y_center(2)), ...
                               gaussmf_val(y_mod(k), fis.sigma, fis.Y_center(3))];
                    for i = 1:fis.rules_number
                        output = output +  degrees(i) * (fis.a_param(i)*sinh(y_mod(k)/36) + fis.b_param(i));
                        w = w + degrees(i);
                    end
                    y_out(k) = output / w;
                end

                % Ograniczenia wartości sygnału wyjściowego, tj. wysokości h_2
                y_out(k) = min(y_out(k), obj.y_max);
                y_out(k) = max(y_out(k), -obj.y_max);
                   
                % Przepisanie sterowań do wektora przeszłych sterowań
                delta_up = [delta_uk(1); delta_up(1:end-1)];
                delta_uz = [u(2, k) - u(2, k-1); delta_uz(1:end-1)];
            
                U_k_1 = ones(obj.Nu,1)*u(1,k-1);
                Y = ones(obj.N,1)*y_out(k);
                Y_0 = Y + obj.M_p*delta_up + obj.M_zP*delta_uz;
                Y_zad = ones(obj.N,1)*y_zad(k);
                B = [-U_min + U_k_1;
                    U_max - U_k_1;
                    -Y_min + Y_0;
                    Y_max - Y_0];
                f = -2*obj.M'*(Y_zad - Y_0);
                
                % Obliczenie przyrostu sterowania dla chwili (i+1) w chwili i-tej
                options = optimoptions('quadprog', 'Display', 'off');
                delta_uk = quadprog(H,f,A,B,[],[],delta_U_min,delta_U_max,[],options);
                delta_u(k) = delta_uk(1);
                
                % Obliczenie sterowania
                u(1,k) = u(1,k-1) + delta_uk(1);
            end
            E_u = obj.lambda .* sum(delta_u.^2)/kk;
            E_y = sum((y_zad - y_out).^2)/kk;
            E =  E_u + E_y;
        end

        function [y_out, u, E, E_u, E_y] = dmc_slWiener(obj, y_zad, u_D, fis, delay, kk, F_10, F_D0, obiekt, type)
            %% Alokacja pamięci
            u = zeros(2, kk);
            u(2,:) = u_D;

            delta_up = zeros(1, obj.D-1);
            delta_uz = zeros(1, obj.D_disturbance);
            delta_uk = 0;
            delta_u = zeros(1, kk);
            
            %% Sterowanie DMC
            y_mod = zeros(size(y_zad));
            y_out = zeros(size(y_zad));
            y_mod(1:delay+2) = y_zad(1:delay+2);
            for k = delay+3:kk
                obiekt.linearization(F_10 + u(1,k-(delay+1)), F_D0 + u(2,k));
                [a, b, s, s_disturbance] = obiekt.mse();
                obj.dynamic_matrix(s');
                obj.past_matrix(s');
                obj.matrix_disturbance(s_disturbance');

                y_mod(k) = - a*y_mod(k-1:-1:k-2)' + b*u(1, k-(delay+1))' + b*u(2, k-1)';

                if (strcmp(type, 'linear'))
                    y_out(k) = evalfis(fis, y_mod(k));
                else
                    gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));

                    output = 0;
                    w = 0;
                    degrees = [gaussmf_val(y_mod(k), fis.sigma, fis.Y_center(1)), ...
                               gaussmf_val(y_mod(k), fis.sigma, fis.Y_center(2)), ...
                               gaussmf_val(y_mod(k), fis.sigma, fis.Y_center(3))];
                    for i = 1:fis.rules_number
                        output = output +  degrees(i) * (fis.a_param(i)*sinh(y_mod(k)/36) + fis.b_param(i));
                        w = w + degrees(i);
                    end
                    y_out(k) = output / w;
                end

                % Ograniczenia wartości sygnału wyjściowego, tj. wysokości h_2
                y_out(k) = min(y_out(k), obj.y_max);
                y_out(k) = max(y_out(k), -obj.y_max);

                % Przepisanie sterowań do wektora przeszłych sterowań
                delta_up = [delta_uk, delta_up(1:end-1)];
                delta_uz = [u(2, k) - u(2, k-1) , delta_uz(1:end-1)];

                % Oblicznie uchybu    
                e = y_zad(k) - y_out(k);
            
                % Obliczenie przyrostu sterowania dla chwili (i+1) w chwili i-tej
                delta_uk = obj.ke * e - obj.ku * delta_up' - obj.kz * delta_uz';
                
                % Ograniczenie wartości przyrostu sterowania
                delta_uk = min(delta_uk, obj.delta_uk_max);
                delta_uk = max(delta_uk, -obj.delta_uk_max);
                delta_u(k) = delta_uk;
                
                % Obliczenie sterowania dla chwili (i+1) w chwili i-tej
                u(1, k) = u(1, k-1) + delta_uk; 
                
                % Ograniczenie sterowania
                u(1, k) = min(u(1, k), obj.u_max);
                u(1, k) = max(u(1, k), -obj.u_max);
            end
            E_u = obj.lambda .* sum(delta_u.^2)/kk;
            E_y = sum((y_zad - y_out).^2)/kk;
            E =  E_u + E_y;
        end

        function [y_mod, u, E, E_u, E_y] = dmc_nplWiener(obj, y_zad, u_D, fis, delay, kk, F_10, F_D0, obiekt, type)
            %% Alokacja pamięci
            u = zeros(2, kk);
            u(2,:) = u_D;

            % delta_up = zeros(1, obj.D-1);
            % delta_uz = zeros(1, obj.D_disturbance);
            % delta_uk = 0;
            delta_u = zeros(1, kk);
            
            %% Sterowanie DMC
            y_mod = zeros(size(y_zad));
            y_out = zeros(size(y_zad));
            y_mod(1:delay+1) = y_zad(1:delay+1);

            h_20 = obiekt.h_20;
            V_1 = zeros(size(y_zad));
            V_2 = zeros(size(y_zad));
            V_1(1:delay+1) = obiekt.V_10;
            V_2(1:delay+1) = obiekt.V_20;
            V_1e = zeros(size(y_zad));
            V_2e = zeros(size(y_zad));
            V_1e(1:delay+1) = obiekt.V_10;
            V_2e(1:delay+1) = obiekt.V_20;

            for k = delay+2:kk

                % Symulacja nieliniowego modelu procesu
                [~, V_1(k), V_2(k), V_1e(k), V_2e(k)] = obiekt.realSimulation(...
                    [u(1, k-(delay+1)) u(2, k-1)]', F_10, F_D0, h_20, ...
                    V_1(k-1), V_2(k-1), V_1e(k-1), V_2e(k-1));

                % Linearyzacja
                obiekt.linearization(F_10 + u(1, k-(delay+1)), F_D0 + u(2, k-1));
                [~, ~, s, ~] = obiekt.mse();

                obj.dynamic_matrix(s');

                % Predykcja swobodnej odpowiedzi
                Y_0 = obiekt.freeAnswer([u(1, k-(delay+1)) u(2, k-1)]', F_10, F_D0, h_20, ...
                                      V_1(k), V_2(k), V_1e(k), V_2e(k), obj.N);

                Y_zad = y_zad(k) * ones(obj.N, 1);
                delta_uk = obj.K * (Y_zad - Y_0);

                % Ograniczenia przyrostów sterowania
                for i = 1:obj.Nu
                    delta_uk = min(delta_uk, obj.delta_uk_max);
                    delta_uk = max(delta_uk, -obj.delta_uk_max);
                end

                % Aktualizacja sterowania
                u(1, k) = u(1, k-1) + delta_uk(1);

                % Ograniczenia sterowania
                u(1, k) = min(u(1, k), obj.u_max);
                u(1, k) = max(u(1, k), -obj.u_max);
                
                % Zapamiętanie przyrostu sterowania
                delta_u(k) = delta_uk(1);
                
                % Predykcja wyjścia
                Y = Y_0 + obj.M * delta_uk;

                y_mod(k) = Y(1);

                if (strcmp(type, 'linear'))
                    y_out(k) = evalfis(fis, y_mod(k));
                else
                    gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));

                    output = 0;
                    w = 0;
                    degrees = [gaussmf_val(y_mod(k), fis.sigma, fis.Y_center(1)), ...
                               gaussmf_val(y_mod(k), fis.sigma, fis.Y_center(2)), ...
                               gaussmf_val(y_mod(k), fis.sigma, fis.Y_center(3))];
                    for i = 1:fis.rules_number
                        output = output +  degrees(i) * (fis.a_param(i)*sinh(y_mod(k)/36) + fis.b_param(i));
                        w = w + degrees(i);
                    end
                    y_out(k) = output / w;
                end

                % Ograniczenia wartości sygnału wyjściowego, tj. wysokości h_2
                y_out(k) = min(y_out(k), obj.y_max);
                y_out(k) = max(y_out(k), -obj.y_max);
            end
            E_u = obj.lambda .* sum(delta_u.^2)/kk;
            E_y = sum((y_zad - y_out).^2)/kk;
            E =  E_u + E_y;
        end

        function [y_out, u, E, E_u, E_y] = dmc_fuzzyWiener(obj, y_zad, u_D, a, b, fis, delay, kk, F_10, F_D0, h_20, obiekt, type)
            gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));
            %% Alokacja pamięci
            u = zeros(2, kk);
            u(2,:) = u_D;

            delta_up = zeros(1, obj.D-1);
            delta_uz = zeros(1, obj.D_disturbance);
            delta_uk = 0;
            delta_u = zeros(1, kk);

            F_1 = [-30 -5 18 30];
            F_D = [-10 -2 10 10];
            L = [1 1 4 4];
            sigma = 12;
            rules_number = 4;
            
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

            disp(Y);
            
            %% Sterowanie DMC
            y_mod = zeros(size(y_zad));
            y_out = zeros(size(y_zad));
            y_mod(1:delay+1) = y_zad(1:delay+1);
            for k = delay+2:kk
                y_mod(k) = - a*y_mod(k-1:-1:k-2)' + b*u(1, k-(delay+1))' + b*u(2, k-1)';

                if (strcmp(type, 'linear'))
                    y_out(k) = evalfis(fis, y_mod(k));
                else
                    gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));

                    output = 0;
                    w = 0;
                    degrees = [gaussmf_val(y_mod(k), fis.sigma, fis.Y_center(1)), ...
                               gaussmf_val(y_mod(k), fis.sigma, fis.Y_center(2)), ...
                               gaussmf_val(y_mod(k), fis.sigma, fis.Y_center(3))];
                    for i = 1:fis.rules_number
                        output = output +  degrees(i) * (fis.a_param(i)*sinh(y_mod(k)/36) + fis.b_param(i));
                        w = w + degrees(i);
                    end
                    y_out(k) = output / w;
                end

                % Ograniczenia wartości sygnału wyjściowego, tj. wysokości h_2
                y_out(k) = min(y_out(k), obj.y_max);
                y_out(k) = max(y_out(k), -obj.y_max);

                % Przepisanie sterowań do wektora przeszłych sterowań
                delta_up = [delta_uk, delta_up(1:end-1)];
                delta_uz = [u(2, k) - u(2, k-1) , delta_uz(1:end-1)];

                % Oblicznie uchybu    
                e = y_zad(k) - y_out(k);

                output = 0;
                w = 0;
                for i = 1:rules_number
                    degree = gaussmf_val(y_out(k), sigma, Y(i));
                    output = output + degree * (KE(i) * e - KU{i} * delta_up' - KZ{i} * delta_uz');
                    w = w + degree;
                end
                delta_uk = output / w;

                % Ograniczenie wartości przyrostu sterowania
                delta_uk = min(delta_uk, obj.delta_uk_max);
                delta_uk = max(delta_uk, -obj.delta_uk_max);
                delta_u(k) = delta_uk;
                
                % Obliczenie sterowania dla chwili (i+1) w chwili i-tej
                u(1, k) = u(1, k-1) + delta_uk;
                
                % Ograniczenie sterowania
                u(1, k) = min(u(1, k), obj.u_max);
                u(1, k) = max(u(1, k), -obj.u_max);
            end
            E_u = obj.lambda .* sum(delta_u.^2)/kk;
            E_y = sum((y_zad - y_mod).^2)/kk;
            E =  E_u + E_y;
        end
        
        function show_result(~, y_mod, y_zad, u, E, E_u, E_y, kk, s, t)
            fprintf("Błąd DMC-%s (%s) \t E_u: %.3f \n", s, t, E_u);
            fprintf("Błąd DMC-%s (%s) \t E_y: %.3f \n", s, t, E_y);
            fprintf("Błąd DMC-%s (%s) \t E = %.3f \n\n", s, t, E);

            figure;
            hold on;
            stairs(0:kk-1, y_mod, 'b-', 'LineWidth', 0.8);
            stairs(0:kk-1, y_zad, 'r-', 'LineWidth', 0.8);
            hold off;
            legend('y', 'y_{zad}', 'Location', 'best');
            title(sprintf('Sygnał wyjściowy y(k) \t DMC - %s \t %s', s, t));
            xlabel('k');
            ylabel('y(k)');
            grid on;
            % saveas(gcf, sprintf('D:/EiTI/MGR/raporty/raport_MGR/pictures/y_%s%s.png', s, t));  % Zapisuje jako plik PNG

            figure;
            hold on;
            stairs(0:kk-1, u(1,:), 'b-', 'LineWidth', 0.8);
            hold off;
            legend('u', 'Location', 'best');
            title(sprintf('Sygnał sterujący u(k) \t DMC - %s \t %s', s, t));
            xlabel('k');
            ylabel('u(k)');
            grid on;
            % saveas(gcf, sprintf('D:/EiTI/MGR/raporty/raport_MGR/pictures/u_%s%s.png', s, t));  % Zapisuje jako plik PNG
        end

    end
end