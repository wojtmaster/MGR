classdef DMC
    properties
        N = 105
        Nu = 2
        D = 105
        M
        M_p
        K
        ke
        ku
        lambda = 0.1
        D_disturbance = 100
        M_zP
        kz
    end

    methods
        function obj = DMC(s, s_disturbance)
            obj = obj.matrix(s);
            obj = obj.matrix_disturbance(s_disturbance);
        end

        function obj = matrix(obj, s)
            obj.M = zeros(obj.N, obj.Nu);
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
            obj.ku = obj.K(1, :) * obj.M_p;
        end

        function obj = matrix_disturbance(obj, s)
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
    
        function dmc_analitic(obj, y_zad, a, b, delay, kk)
            %% Alokacja pamięci
            u = zeros(2, kk);
            u(2,:) = [repelem((rand(1, kk/100) * 10 - 5), 100)];

            delta_up = zeros(1, obj.D-1);
            delta_uz = zeros(1, obj.D_disturbance);
            delta_uk = 0;
            
            %% Ograniczenia sygnału sterującego i wyjścia obiektu
            u_max = 50;
            delta_uk_max = 10;
            y_max = 30;
            
            %% Sterowanie DMC
            y_mod = zeros(size(y_zad));
            y_mod(1:delay+2) = y_zad(1:delay+2);
            for k = delay+3:kk
                y_mod(k) = - a*[y_mod(k-1:-1:k-2)]' + b*[u(1, k-(delay+1):-1:k-(delay+2))]' + b*[u(2, k-1:-1:k-2)]';

                % Ograniczenia wartości sygnału wyjściowego, tj. wysokości h_2
                y_mod(k) = min(y_mod(k), y_max);
                y_mod(k) = max(y_mod(k), -y_max);

                % Przepisanie sterowań do wektora przeszłych sterowań
                delta_up = [delta_uk, delta_up(1:end-1)];
                delta_uz = [u(2, k) - u(2, k-1) , delta_uz(1:end-1)];

                % Oblicznie uchybu    
                e = y_zad(k) - y_mod(k);
            
                % Obliczenie przyrostu sterowania dla chwili (i+1) w chwili i-tej
                % delta_uk = obj.ke * e - obj.ku * delta_up' - obj.kz * delta_uz';
                delta_uk = obj.ke * e - obj.ku * delta_up';
                
                % Ograniczenie wartości przyrostu sterowania
                delta_uk = min(delta_uk, delta_uk_max);
                delta_uk = max(delta_uk, -delta_uk_max);
                
                % Obliczenie sterowania dla chwili (i+1) w chwili i-tej
                u(1, k) = u(1, k-1) + delta_uk; 
                
                % Ograniczenie sterowania
                u(1, k) = min(u(1, k), u_max);
                u(1, k) = max(u(1, k), -u_max);
            end

            E = sum((y_zad - y_mod).^2)/(kk);
            fprintf("Błąd DMC \n");
            fprintf("E = %.3f \n", E);
            
            figure;
            hold on;
            stairs(0:kk-1, y_mod, 'b-', 'LineWidth', 0.8);
            stairs(0:kk-1, y_zad, 'r-', 'LineWidth', 0.8);
            hold off;
            legend('y_{mod}', 'y_{zad}');
            title('y(k)');
            xlabel('k');
            ylabel('y');
            grid on;
            
            figure;
            hold on;
            stairs(0:kk-1, u(1,:), 'b-', 'LineWidth', 0.8);
            hold off;
            legend('u');
            title('u(k)');
            xlabel('k');
            ylabel('u');
            grid on;

        end
    
        function dmc_numeric(obj, y_zad, a, b, delay, kk)
            % Alokacja pamięci
            u = zeros(2, kk);
            u(2,:) = [repelem((rand(1, kk/100) * 10 - 5), 100)];
            
            delta_up = zeros(obj.D-1,1);
            delta_uk = zeros(obj.Nu,1);
            delta_uz = zeros(obj.D_disturbance, 1);

            
            u_max = 50;
            delta_uk_max = 10;
            y_max = 30;

            Y_max = ones(obj.N,1)*y_max;
            Y_min = ones(obj.N,1)*(-y_max);
            
            U_max = ones(obj.Nu,1)*u_max;
            U_min = ones(obj.Nu,1)*(-u_max);
            delta_U_max = ones(obj.Nu,1)*delta_uk_max;
            delta_U_min = ones(obj.Nu,1)*(-delta_uk_max);
            
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
                % Y_0 = Y + obj.M_p*delta_up + obj.M_zP*delta_uz;
                Y_0 = Y + obj.M_p*delta_up;
                Y_zad = ones(obj.N,1)*y_zad(k);
                B = [-U_min + U_k_1;
                    U_max - U_k_1;
                    -Y_min + Y_0;
                    Y_max - Y_0];
                f = -2*obj.M'*(Y_zad - Y_0);
                
                % Obliczenie przyrostu sterowania dla chwili (i+1) w chwili i-tej
                delta_uk = quadprog(H,f,A,B,[],[],delta_U_min,delta_U_max);
                
                % Obliczenie sterowania
                u(1,k) = u(1,k-1) + delta_uk(1); 
            end

            E = sum((y_zad - y_mod).^2)/(kk);
            fprintf("Błąd DMC \n");
            fprintf("E = %.3f \n", E);
            
            figure;
            hold on;
            stairs(0:kk-1, y_mod, 'b-', 'LineWidth', 0.8);
            stairs(0:kk-1, y_zad, 'r-', 'LineWidth', 0.8);
            hold off;
            legend('y_{mod}', 'y_{zad}');
            title('y(k)');
            xlabel('k');
            ylabel('y');
            grid on;
            
            figure;
            hold on;
            stairs(0:kk-1, u(1,:), 'b-', 'LineWidth', 0.8);
            hold off;
            legend('u');
            title('u(k)');
            xlabel('k');
            ylabel('u');
            grid on;
        end
    
        function dmc_sl(obj)
            disp(obj.N);
        end
    end
end

























