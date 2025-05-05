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
    
        function [y_mod, u, E, E_u, E_y] = dmc_analitic(obj, y_zad, u_D, a, b, delay, kk)
            %% Alokacja pamięci
            u = zeros(2, kk);
            u(2,:) = u_D;

            delta_up = zeros(1, obj.D-1);
            delta_uz = zeros(1, obj.D_disturbance);
            delta_uk = 0;
            delta_u = zeros(1, kk);
            
            %% Sterowanie DMC
            y_mod = zeros(size(y_zad));
            y_mod(1:delay+2) = y_zad(1:delay+2);
            for k = delay+3:kk
                y_mod(k) = - a*[y_mod(k-1:-1:k-2)]' + b*[u(1, k-(delay+1):-1:k-(delay+2))]' + b*[u(2, k-1:-1:k-2)]';

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
    
        function [y_mod, u, E] = dmc_numeric(obj, y_zad, u_D, a, b, delay, kk)
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
            y_mod(1:delay+2) = y_zad(1:delay+2);
            for k = delay+3:kk
                y_mod(k) = - a*[y_mod(k-1:-1:k-2)]' + b*[u(1, k-(delay+1):-1:k-(delay+2))]' + b*[u(2, k-1:-1:k-2)]';
                   
                % Przepisanie sterowań do wektora przeszłych sterowań
                delta_up = [delta_uk(1); delta_up(1:end-1)];
                delta_uz = [u(2, k) - u(2, k-1); delta_uz(1:end-1)];
            
                U_k_1 = ones(obj.Nu,1)*u(1,k-1);
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
                delta_u(k) = delta_uk(1);
                
                % Obliczenie sterowania
                u(1,k) = u(1,k-1) + delta_uk(1);
            end
            E = sum((y_zad - y_mod).^2)/kk + obj.lambda .* (sum(delta_u.^2))/kk;
        end
    
        function [y_mod, u, E] = dmc_sl(obj, y_zad, u_D, delay, kk, F_10, F_D0, linearization)
            %% Alokacja pamięci
            u = zeros(2, kk);
            u(2,:) = u_D;

            delta_up = zeros(1, obj.D-1);
            delta_uz = zeros(1, obj.D_disturbance);
            delta_uk = 0;
            delta_u = zeros(1, kk);
            
            %% Sterowanie DMC
            y_mod = zeros(size(y_zad));
            y_mod(1:delay+2) = y_zad(1:delay+2);
            for k = delay+3:kk
                [s, s_disturbance, a, b] = linearization(F_10 + u(1,k), F_D0 + u(2,k));
                obj.dynamic_matrix(s);
                obj.past_matrix(s);
                obj.matrix_disturbance(s_disturbance);

                y_mod(k) = - a*[y_mod(k-1:-1:k-2)]' + b*[u(1, k-(delay+1):-1:k-(delay+2))]' + b*[u(2, k-1:-1:k-2)]';

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
                delta_u(k) = delta_uk;
                
                % Obliczenie sterowania dla chwili (i+1) w chwili i-tej
                u(1, k) = u(1, k-1) + delta_uk; 
                
                % Ograniczenie sterowania
                u(1, k) = min(u(1, k), obj.u_max);
                u(1, k) = max(u(1, k), -obj.u_max);
            end
            E = sum((y_zad - y_mod).^2)/kk + obj.lambda .* (sum(delta_u.^2))/kk;
        end

        function [y_mod, u, E] = dmc_SL(obj, y_zad, u_D, delay, kk, F_10, F_D0, linearization)
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
            
            y_mod = zeros(size(y_zad));
            y_mod(1:delay+2) = y_zad(1:delay+2);
            for k = delay+3:kk
                [s, s_disturbance, a, b] = linearization(F_10 + u(1,k), F_D0 + u(2,k));
                obj.dynamic_matrix(s);
                obj.past_matrix(s);
                obj.matrix_disturbance(s_disturbance);

                y_mod(k) = - a*[y_mod(k-1:-1:k-2)]' + b*[u(1, k-(delay+1):-1:k-(delay+2))]' + b*[u(2, k-1:-1:k-2)]';
                   
                % Przepisanie sterowań do wektora przeszłych sterowań
                delta_up = [delta_uk(1); delta_up(1:end-1)];
                delta_uz = [u(2, k) - u(2, k-1); delta_uz(1:end-1)];
            
                A = [-J; J; -obj.M; obj.M];
                H = 2*(obj.M'*obj.M + obj.lambda*eye(obj.Nu));
                U_k_1 = ones(obj.Nu,1)*u(1,k-1);
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
                delta_u(k) = delta_uk(1);
                
                % Obliczenie sterowania
                u(1,k) = u(1,k-1) + delta_uk(1); 
            end
            E = sum((y_zad - y_mod).^2)/kk + obj.lambda .* (sum(delta_u.^2))/kk;
        end

        function [y_mod, u, E] = dmc_npl(obj, y_zad, u_D, delay, kk, F_10, F_D0, linearization, rk4)
            %% Alokacja pamięci
            u = zeros(2, kk);
            u(2,:) = u_D;

            delta_up = zeros(1, obj.D-1);
            delta_uz = zeros(1, obj.D_disturbance);
            delta_uk = 0;
            delta_u = zeros(1, kk);
            
            %% Sterowanie DMC
            y_mod = zeros(size(y_zad));
            y_mod(1:delay+2) = y_zad(1:delay+2);
            for k = delay+3:kk
                [s, ~, a, b] = linearization(F_10 + u(1, k), F_D0 + u(2, k));
                [S, ~] = rk4([ones(1, obj.D); zeros(1, obj.D)], obj.D);
                [S_disturbance, ~] = rk4([zeros(1, obj.D_disturbance); ones(1, obj.D_disturbance)], obj.D_disturbance);

                obj.dynamic_matrix(s);
                obj.past_matrix(S');
                obj.matrix_disturbance(S_disturbance');

                y_mod(k) = - a*[y_mod(k-1:-1:k-2)]' + b*[u(1, k-(delay+1):-1:k-(delay+2))]' + b*[u(2, k-1:-1:k-2)]';

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
                delta_u(k) = delta_uk;
                
                % Obliczenie sterowania dla chwili (i+1) w chwili i-tej
                u(1, k) = u(1, k-1) + delta_uk; 
                
                % Ograniczenie sterowania
                u(1, k) = min(u(1, k), obj.u_max);
                u(1, k) = max(u(1, k), -obj.u_max);
            end
            E = sum((y_zad - y_mod).^2)/kk + obj.lambda .* (sum(delta_u.^2))/kk;
        end

        function [y_mod, u, E] = dmc_NPL(obj, y_zad, u_D, delay, kk, F_10, F_D0, linearization, rk4)
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
            
            y_mod = zeros(size(y_zad));
            y_mod(1:delay+2) = y_zad(1:delay+2);
            for k = delay+3:kk
                [s, ~, a, b] = linearization(F_10 + u(1, k), F_D0 + u(2, k));
                [S, ~] = rk4([ones(1, obj.D); zeros(1, obj.D)], obj.D);
                [S_disturbance, ~] = rk4([zeros(1, obj.D_disturbance); ones(1, obj.D_disturbance)], obj.D_disturbance);

                obj.dynamic_matrix(s);
                obj.past_matrix(S');
                obj.matrix_disturbance(S_disturbance');

                y_mod(k) = - a*[y_mod(k-1:-1:k-2)]' + b*[u(1, k-(delay+1):-1:k-(delay+2))]' + b*[u(2, k-1:-1:k-2)]';
                   
                % Przepisanie sterowań do wektora przeszłych sterowań
                delta_up = [delta_uk(1); delta_up(1:end-1)];
                delta_uz = [u(2, k) - u(2, k-1); delta_uz(1:end-1)];
            
                A = [-J; J; -obj.M; obj.M];
                H = 2*(obj.M'*obj.M + obj.lambda*eye(obj.Nu));
                U_k_1 = ones(obj.Nu,1)*u(1,k-1);
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
                delta_u(k) = delta_uk(1);
                
                % Obliczenie sterowania
                u(1,k) = u(1,k-1) + delta_uk(1); 
            end
            E = sum((y_zad - y_mod).^2)/kk + obj.lambda .* (sum(delta_u.^2))/kk;
        end

        function [y_mod, u, E] = dmc_fuzzy(obj, y_zad, u_D, a, b, delay, kk, fis, F_10, F_D0, linearization)
            %% Alokacja pamięci
            u = zeros(2, kk);
            u(2,:) = u_D;

            delta_up = zeros(1, obj.D-1);
            delta_uz = zeros(1, obj.D_disturbance);
            delta_uk = 0;
            delta_u = zeros(1, kk);

            F_1 = [-45, 45];
            F_D = [-5, 5];
            
            KE = zeros(1, length(F_1)*length(F_D));
            KU = cell(1, length(F_1)*length(F_D));
            KZ = cell(1, length(F_1)*length(F_D));

            for i = 1:length(F_1)
                for j = 1:length(F_D)
                    [s, s_disturbance, ~, ~] = linearization(F_10+F_1(i), F_D0+F_D(j));
                    obj.dynamic_matrix(s);
                    obj.past_matrix(s);
                    obj.matrix_disturbance(s_disturbance);
            
                    KE((i - 1) * length(F_D) + j) = obj.ke;
                    KU{(i - 1) * length(F_D) + j} = obj.ku;
                    KZ{(i - 1) * length(F_D) + j} = obj.kz;
                end
            end
            
            %% Sterowanie DMC
            y_mod = zeros(size(y_zad));
            y_mod(1:delay+2) = y_zad(1:delay+2);
            for k = delay+3:kk
                y_mod(k) = - a*[y_mod(k-1:-1:k-2)]' + b*[u(1, k-(delay+1):-1:k-(delay+2))]' + b*[u(2, k-1:-1:k-2)]';

                % Ograniczenia wartości sygnału wyjściowego, tj. wysokości h_2
                y_mod(k) = min(y_mod(k), obj.y_max);
                y_mod(k) = max(y_mod(k), -obj.y_max);

                % Przepisanie sterowań do wektora przeszłych sterowań
                delta_up = [delta_uk, delta_up(1:end-1)];
                delta_uz = [u(2, k) - u(2, k-1) , delta_uz(1:end-1)];

                % Oblicznie uchybu    
                e = y_zad(k) - y_mod(k);

                % Obliczenie przyrostu sterowania dla chwili (i+1) w chwili i-tej
                inputs = [y_mod(k-1), u(1, k-(delay+1)), u(2, k-1)];
                % inputs = [y_mod(k-1), u(1, k-(delay+1))];
                [~, degree] = evalfis(fis, inputs);
                strength = max(degree);
                % strength = strength / sum(strength);
                for i = 1:length(strength)
                    idx = mod(i-1, length(F_1)*length(F_D)) + 1;
                    delta_uk = delta_uk + strength(i)*(KE(idx) * e - KU{idx} * delta_up' - KZ{idx} * delta_uz');
                end
                delta_uk = delta_uk / sum(strength);

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
            E = sum((y_zad - y_mod).^2)/kk + obj.lambda .* sum(delta_u.^2)/kk;
        end

        function [y_mod, u, E] = dmc_FUZZY(obj, y_zad, u_D, a, b, delay, kk, fis, F_10, F_D0, linearization)
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

            F_1 = [-45, 45];
            F_D = [-5, 5];
            
            m = cell(1, length(F_1)*length(F_D));
            m_p = cell(1, length(F_1)*length(F_D));
            m_zP = cell(1, length(F_1)*length(F_D));
            
            for i = 1:length(F_1)
                for j = 1:length(F_D)
                    [s, s_disturbance, ~, ~] = linearization(F_10+F_1(i), F_D0+F_D(j));
                    obj.dynamic_matrix(s);
                    obj.past_matrix(s);
                    obj.matrix_disturbance(s_disturbance);
            
                    m{(i - 1) * length(F_D) + j} = obj.M;
                    m_p{(i - 1) * length(F_D) + j} = obj.M_p;
                    m_zP{(i - 1) * length(F_D) + j} = obj.M_zP;
                end
            end
            
            y_mod = zeros(size(y_zad));
            y_mod(1:delay+2) = y_zad(1:delay+2);
            for k = delay+3:kk
                y_mod(k) = - a*[y_mod(k-1:-1:k-2)]' + b*[u(1, k-(delay+1):-1:k-(delay+2))]' + b*[u(2, k-1:-1:k-2)]';
                   
                % Przepisanie sterowań do wektora przeszłych sterowań
                delta_up = [delta_uk(1); delta_up(1:end-1)];
                delta_uz = [u(2, k) - u(2, k-1); delta_uz(1:end-1)];
            
                U_k_1 = ones(obj.Nu,1)*u(1,k-1);
                Y = ones(obj.N,1)*y_mod(k);
                Y_zad = ones(obj.N,1)*y_zad(k);

                inputs = [y_mod(k-1), u(1, k-(delay+1)), u(2, k-1)];
                [~, degree] = evalfis(fis, inputs);
                % strength = prod(degree,2);
                strength = max(degree);
                for i = 1:length(strength)
                    idx = mod(i-1, length(F_1)*length(F_D)) + 1;

                    A = [-J; J; -m{idx}; m{idx}];
                    H = 2*(m{idx}'*m{idx} + obj.lambda*eye(obj.Nu));

                    Y_0 = Y + m_p{idx}*delta_up + m_zP{idx}*delta_uz;
                    f = -2*m{idx}'*(Y_zad - Y_0);

                    B = [-U_min + U_k_1;
                    U_max - U_k_1;
                    -Y_min + Y_0;
                    Y_max - Y_0];

                    % Obliczenie przyrostu sterowania dla chwili (i+1) w chwili i-tej
                    options = optimoptions('quadprog', 'Display', 'off');
                    delta_uk = delta_uk + strength(i).*quadprog(H,f,A,B,[],[],delta_U_min,delta_U_max,[],options);
                end
                delta_uk = delta_uk / sum(strength);
                
                % Obliczenie sterowania
                u(1,k) = u(1,k-1) + delta_uk(1);
            end
            E = sum((y_zad - y_mod).^2)/kk + obj.lambda .* (sum(delta_u.^2))/kk;
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
                    y_fuzzy(k-(delay+1)) = evalfis(fis, [u(1,k-(delay+1)), u(2,k-1)]);
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
                    y_fuzzy(k-(delay+1)) = output / w;
                end
                
                if(u(1,k-(delay+1)) ~= 0)
                    gain = y_fuzzy(k-(delay+1)) / (u(1,k-(delay+1)) + u(2,k-1));
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
                    y_fuzzy(k-(delay+1)) = evalfis(fis, [u(1,k-(delay+1)), u(2,k-1)]);
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
                    y_fuzzy(k-(delay+1)) = output / w;
                end

                if(u(1,k-(delay+1)) ~= 0)
                    gain = y_fuzzy(k-(delay+1)) / (u(1,k-(delay+1)) + u(2,k-1));
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

        function [y_mod, u, E, E_u, E_y] = dmc_slHammerstein(obj, y_zad, u_D, a, b, fis, delay, kk, F_10, F_D0, h_20, linearization, obiekt, type)
            %% Alokacja pamięci
            u = zeros(2, kk);
            u(2,:) = u_D;

            delta_up = zeros(1, obj.D-1);
            delta_uz = zeros(1, obj.D_disturbance);
            delta_uk = 0;
            delta_u = zeros(1, kk);
            
            %% Sterowanie DMC
            y_mod = zeros(size(y_zad));
            y_mod(1:delay+2) = y_zad(1:delay+2);
            y_fuzzy = zeros(size(y_zad));
            for k = delay+3:kk
                % [s, s_disturbance] = obiekt.linearization(F_10 + u(1,k-1), F_D0 + u(2,k-1));
                % obj.dynamic_matrix(s);
                % obj.past_matrix(s);
                % obj.matrix_disturbance(s_disturbance);

                delta_F10 = obiekt.F_10 - F_10;
                delta_FD0 = obiekt.F_D0 - F_D0;
                delta_h20 = obiekt.h_20 - h_20;
                
                if (mod(k,200) == 0)
                    [s, s_disturbance] = obiekt.linearization(F_10 + u(1,k-1), F_D0 + u(2,k-1));
                    h_20 = obiekt.h_20;

                    obj.dynamic_matrix(s);
                    obj.past_matrix(s);
                    obj.matrix_disturbance(s_disturbance);
                    [a, b] = obiekt.sopdt();
                end

                if (strcmp(type, 'linear'))
                    y_fuzzy(k-(delay+2)) = evalfis(fis, [u(1,k-(delay+2)), u(2,k-2)]);
                else
                    gaussmf_val = @(x, sigma, c) exp(-((x - c).^2) / (2 * sigma^2));

                    deg_u1 = [gaussmf_val(u(1,k-(delay+2)), fis.sigma_F1, fis.F1_center(1)), gaussmf_val(u(1,k-(delay+2)), fis.sigma_F1, fis.F1_center(2)), gaussmf_val(u(1,k-(delay+2)), fis.sigma_F1, fis.F1_center(3))];
                    deg_u2 = [gaussmf_val(u(2,k-2), fis.sigma_FD, fis.FD_center(1)), gaussmf_val(u(2,k-2), fis.sigma_FD, fis.FD_center(2))];
                    
                    degrees_all(1) = deg_u1(1) * deg_u2(1);
                    degrees_all(2) = deg_u1(1) * deg_u2(2);
                    degrees_all(3) = deg_u1(2) * deg_u2(1);
                    degrees_all(4) = deg_u1(2) * deg_u2(2);
                    degrees_all(5) = deg_u1(3) * deg_u2(1);
                    degrees_all(6) = deg_u1(3) * deg_u2(2);
                    
                    w = 0;
                    output = 0;
                    for i = 1:fis.rules_number
                        output = output + degrees_all(i)*(fis.a_param(i)*sinh(u(1,k-(delay+2))/22.5) + fis.b_param(i)*sinh(u(2,k-2)/7.5) + fis.c_param(i));
                        w = w + degrees_all(i);
                    end
                    y_fuzzy(k-(delay+2)) = output / w;
                end

                if(u(1,k-(delay+2)) ~= 0)
                    gain = (y_fuzzy(k-(delay+2)) - h_20) / (u(1,k-(delay+2)) + u(2,k-2));
                else
                    gain = 1;
                end

                y_mod(k) = - a*[y_mod(k-1:-1:k-2)]' + gain*b*[u(1, k-(delay+1):-1:k-(delay+2))]' + gain*b*[u(2, k-1:-1:k-2)]';

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

        function [y_mod, u_fuzzy, E, E_u, E_y] = dmc_nplHammerstein(obj, y_zad, u_D, fis, delay, kk, F_10, F_D0, linearization, rk4, type)
            %% Alokacja pamięci
            u = zeros(2, kk);
            u_fuzzy = zeros(1, kk);
            u(2,:) = u_D;

            delta_up = zeros(1, obj.D-1);
            delta_uz = zeros(1, obj.D_disturbance);
            delta_uk = 0;
            delta_u = zeros(1, kk);
            
            %% Sterowanie DMC
            y_mod = zeros(size(y_zad));
            y_mod(1:delay+2) = y_zad(1:delay+2);
            for k = delay+3:kk
                [s, ~, a, b] = linearization(F_10 + u(1, k), F_D0 + u(2, k));
                [S, ~] = rk4([ones(1, obj.D); zeros(1, obj.D)], obj.D);
                [S_disturbance, ~] = rk4([zeros(1, obj.D_disturbance); ones(1, obj.D_disturbance)], obj.D_disturbance);

                obj.dynamic_matrix(s);
                obj.past_matrix(S');
                obj.matrix_disturbance(S_disturbance');

                y_mod(k) = - a*[y_mod(k-1:-1:k-2)]' + b*[u_fuzzy(k-(delay+1):-1:k-(delay+2))]' + b*[u(2, k-1:-1:k-2)]';

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

                if (strcmp(type, 'linear'))
                    u_fuzzy(k) = evalfis(fis, u(1, k));
                else
                    [~, degrees] = evalfis(fis, u(1, k));
                    output = 0;
                    for i = 1:length(fis.Rules)
                        output = output + degrees(i) * (fis.Outputs.MembershipFunctions(i).Parameters(1)*sinh(u(1, k)/fis.Outputs.MembershipFunctions(i).Parameters(2)));
                    end
                    u_fuzzy(k) = output / sum(degrees);
                end
                delta_u(k) = u_fuzzy(k) - u_fuzzy(k-1);
            end
            E_u = obj.lambda .* sum(delta_u.^2)/kk;
            E_y = sum((y_zad - y_mod).^2)/kk;
            E =  E_u + E_y;
        end

        function [y_mod, u_fuzzy, E, E_u, E_y] = dmc_fuzzyHammerstein(obj, y_zad, u_D, a, b, fis, delay, kk, FIS, F_10, F_D0, linearization, type)
            %% Alokacja pamięci
            u = zeros(2, kk);
            u_fuzzy = zeros(1, kk);
            u(2,:) = u_D;

            delta_up = zeros(1, obj.D-1);
            delta_uz = zeros(1, obj.D_disturbance);
            delta_uk = 0;
            delta_u = zeros(1, kk);

            F_1 = [-45, 45];
            F_D = [-5, 5];
            
            KE = zeros(1, length(F_1)*length(F_D));
            KU = cell(1, length(F_1)*length(F_D));
            KZ = cell(1, length(F_1)*length(F_D));

            for i = 1:length(F_1)
                for j = 1:length(F_D)
                    [s, s_disturbance, ~, ~] = linearization(F_10+F_1(i), F_D0+F_D(j));
                    obj.dynamic_matrix(s);
                    obj.past_matrix(s);
                    obj.matrix_disturbance(s_disturbance);
            
                    KE((i - 1) * length(F_D) + j) = obj.ke;
                    KU{(i - 1) * length(F_D) + j} = obj.ku;
                    KZ{(i - 1) * length(F_D) + j} = obj.kz;
                end
            end
            
            %% Sterowanie DMC
            y_mod = zeros(size(y_zad));
            y_mod(1:delay+2) = y_zad(1:delay+2);
            for k = delay+3:kk
                y_mod(k) = - a*[y_mod(k-1:-1:k-2)]' + b*[u_fuzzy(k-(delay+1):-1:k-(delay+2))]' + b*[u(2, k-1:-1:k-2)]';

                % Ograniczenia wartości sygnału wyjściowego, tj. wysokości h_2
                y_mod(k) = min(y_mod(k), obj.y_max);
                y_mod(k) = max(y_mod(k), -obj.y_max);

                % Przepisanie sterowań do wektora przeszłych sterowań
                delta_up = [delta_uk, delta_up(1:end-1)];
                delta_uz = [u(2, k) - u(2, k-1) , delta_uz(1:end-1)];

                % Oblicznie uchybu    
                e = y_zad(k) - y_mod(k);

                % Obliczenie przyrostu sterowania dla chwili (i+1) w chwili i-tej
                inputs = [y_mod(k-1), u(1, k-(delay+1)), u(2, k-1)];
                % inputs = [y_mod(k-1), u(1, k-(delay+1))];
                [~, degree] = evalfis(FIS, inputs);
                strength = max(degree);
                % strength = strength / sum(strength);
                for i = 1:length(strength)
                    idx = mod(i-1, length(F_1)*length(F_D)) + 1;
                    delta_uk = delta_uk + strength(i)*(KE(idx) * e - KU{idx} * delta_up' - KZ{idx} * delta_uz');
                end
                delta_uk = delta_uk / sum(strength);

                % Ograniczenie wartości przyrostu sterowania
                delta_uk = min(delta_uk, obj.delta_uk_max);
                delta_uk = max(delta_uk, -obj.delta_uk_max);
                
                % Obliczenie sterowania dla chwili (i+1) w chwili i-tej
                u(1, k) = u(1, k-1) + delta_uk;
                
                % Ograniczenie sterowania
                u(1, k) = min(u(1, k), obj.u_max);
                u(1, k) = max(u(1, k), -obj.u_max);

                if (strcmp(type, 'linear'))
                    u_fuzzy(k) = evalfis(fis, u(1, k));
                else
                    [~, degrees] = evalfis(fis, u(1, k));
                    output = 0;
                    for i = 1:length(fis.Rules)
                        output = output + degrees(i) * (fis.Outputs.MembershipFunctions(i).Parameters(1)*sinh(u(1, k)/fis.Outputs.MembershipFunctions(i).Parameters(2)));
                    end
                    u_fuzzy(k) = output / sum(degrees);
                end
                delta_u(k) = u_fuzzy(k) - u_fuzzy(k-1);
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

        function [y_out, u, E, E_u, E_y] = dmc_slWiener(obj, y_zad, u_D, fis, delay, kk, F_10, F_D0, linearization, type)
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
                [s, s_disturbance, a, b] = linearization(F_10 + u(1,k), F_D0 + u(2,k));
                obj.dynamic_matrix(s);
                obj.past_matrix(s);
                obj.matrix_disturbance(s_disturbance);

                y_mod(k) = - a*[y_mod(k-1:-1:k-2)]' + b*[u(1, k-(delay+1):-1:k-(delay+2))]' + b*[u(2, k-1:-1:k-2)]';

                if (strcmp(type, 'linear'))
                    y_out(k) = evalfis(fis, y_mod(k));
                else
                    [~, degrees] = evalfis(fis, y_mod(k)); % Przepuszczenie przez model TS
                    output = 0;
                    for i = 1:length(fis.Rules)
                        output = output + degrees(i) * (fis.Outputs.MembershipFunctions(i).Parameters(1)*sinh(y_mod(k)/fis.Outputs.MembershipFunctions(i).Parameters(2)));
                    end
                    y_out(k) = output / sum(degrees);
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
            E_y = sum((y_zad - y_mod).^2)/kk;
            E =  E_u + E_y;
        end

        function [y_out, u, E, E_u, E_y] = dmc_nplWiener(obj, y_zad, u_D, fis, delay, kk, F_10, F_D0, linearization, rk4, type)
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
                [s, ~, a, b] = linearization(F_10 + u(1, k), F_D0 + u(2, k));
                [S, ~] = rk4([ones(1, obj.D); zeros(1, obj.D)], obj.D);
                [S_disturbance, ~] = rk4([zeros(1, obj.D_disturbance); ones(1, obj.D_disturbance)], obj.D_disturbance);

                obj.dynamic_matrix(s);
                obj.past_matrix(S');
                obj.matrix_disturbance(S_disturbance');

                y_mod(k) = - a*[y_mod(k-1:-1:k-2)]' + b*[u(1, k-(delay+1):-1:k-(delay+2))]' + b*[u(2, k-1:-1:k-2)]';

                if (strcmp(type, 'linear'))
                    y_out(k) = evalfis(fis, y_mod(k));
                else
                    [~, degrees] = evalfis(fis, y_mod(k)); % Przepuszczenie przez model TS
                    output = 0;
                    for i = 1:length(fis.Rules)
                        output = output + degrees(i) * (fis.Outputs.MembershipFunctions(i).Parameters(1)*sinh(y_mod(k)/fis.Outputs.MembershipFunctions(i).Parameters(2)));
                    end
                    y_out(k) = output / sum(degrees);
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
            E_y = sum((y_zad - y_mod).^2)/kk;
            E =  E_u + E_y;
        end

        function [y_out, u, E, E_u, E_y] = dmc_fuzzyWiener(obj, y_zad, u_D, a, b, fis, delay, kk, FIS, F_10, F_D0, linearization, type)
            %% Alokacja pamięci
            u = zeros(2, kk);
            u(2,:) = u_D;

            delta_up = zeros(1, obj.D-1);
            delta_uz = zeros(1, obj.D_disturbance);
            delta_uk = 0;
            delta_u = zeros(1, kk);

            F_1 = [-45, 45];
            F_D = [-5, 5];
            
            KE = zeros(1, length(F_1)*length(F_D));
            KU = cell(1, length(F_1)*length(F_D));
            KZ = cell(1, length(F_1)*length(F_D));

            for i = 1:length(F_1)
                for j = 1:length(F_D)
                    [s, s_disturbance, ~, ~] = linearization(F_10+F_1(i), F_D0+F_D(j));
                    obj.dynamic_matrix(s);
                    obj.past_matrix(s);
                    obj.matrix_disturbance(s_disturbance);
            
                    KE((i - 1) * length(F_D) + j) = obj.ke;
                    KU{(i - 1) * length(F_D) + j} = obj.ku;
                    KZ{(i - 1) * length(F_D) + j} = obj.kz;
                end
            end
            
            %% Sterowanie DMC
            y_mod = zeros(size(y_zad));
            y_out = zeros(size(y_zad));
            y_mod(1:delay+2) = y_zad(1:delay+2);
            for k = delay+3:kk
                y_mod(k) = - a*[y_mod(k-1:-1:k-2)]' + b*[u(1, k-(delay+1):-1:k-(delay+2))]' + b*[u(2, k-1:-1:k-2)]';

                if (strcmp(type, 'linear'))
                    y_out(k) = evalfis(fis, y_mod(k));
                else
                    [~, degrees] = evalfis(fis, y_mod(k)); % Przepuszczenie przez model TS
                    output = 0;
                    for i = 1:length(fis.Rules)
                        output = output + degrees(i) * (fis.Outputs.MembershipFunctions(i).Parameters(1)*sinh(y_mod(k)/fis.Outputs.MembershipFunctions(i).Parameters(2)));
                    end
                    y_out(k) = output / sum(degrees);
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
                inputs = [y_mod(k-1), u(1, k-(delay+1)), u(2, k-1)];
                % inputs = [y_mod(k-1), u(1, k-(delay+1))];
                [~, degree] = evalfis(FIS, inputs);
                strength = max(degree);
                % strength = strength / sum(strength);
                for i = 1:length(strength)
                    idx = mod(i-1, length(F_1)*length(F_D)) + 1;
                    delta_uk = delta_uk + strength(i)*(KE(idx) * e - KU{idx} * delta_up' - KZ{idx} * delta_uz');
                end
                delta_uk = delta_uk / sum(strength);

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