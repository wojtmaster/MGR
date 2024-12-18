classdef Hammerstein < Base
    properties
        u
        u_ucz
        u_wer
        u_start
        u_end
        y
        y_ucz
        y_mod_ucz
        y_wer
        y_mod_wer
        R
        params

        u_disturbance
        u_ucz_disturbance
        u_wer_disturbance
        u_start_disturbance
        u_end_disturbance
        y_disturbance
        y_ucz_disturbance
        y_mod_ucz_disturbance;
        y_wer_disturbance
        y_mod_wer_disturbance;
        R_disturbance
        params_disturbance

        n
        mode
        s
    end

    methods 
        function obj = Hammerstein(F_1_start, F_1_end, F_10, F_D_start, F_D_end, F_D0, h_20, alpha_2, n, mode, s)
            F_1 = linspace(F_1_start, F_1_end, n);
            F_D = F_D0;
            h_2 = ((F_1+F_D)/alpha_2).^2;

            obj.y = h_2 - h_20;
            obj.u = F_1 - F_10;
            obj.u_ucz = obj.u(1:2:end-1);
            obj.u_wer = obj.u(2:2:end);
            obj.u_start = F_1_start - F_10;
            obj.u_end = F_1_end - F_10;
            obj.y_ucz = obj.y(1:2:end-1);
            obj.y_mod_ucz = zeros(size(obj.y_ucz));
            obj.y_wer = obj.y(2:2:end);
            obj.y_mod_wer = zeros(size(obj.y_wer));

            F_1 = F_10;
            F_D = linspace(F_D_start, F_D_end, n);
            h_2 = ((F_1+F_D)/alpha_2).^2;

            obj.y_disturbance = h_2 - h_20;
            obj.u_disturbance = F_D - F_D0;
            obj.u_ucz_disturbance = obj.u_disturbance(1:2:end-1);
            obj.u_wer_disturbance = obj.u_disturbance(2:2:end);
            obj.u_start_disturbance = F_D_start - F_D0;
            obj.u_end_disturbance = F_D_end - F_D0;
            obj.y_ucz_disturbance = obj.y_disturbance(1:2:end-1);
            obj.y_mod_ucz_disturbance = zeros(size(obj.y_ucz_disturbance));
            obj.y_wer_disturbance = obj.y_disturbance(2:2:end);
            obj.y_mod_wer_disturbance = zeros(size(obj.y_wer_disturbance));

            obj.n = n;
            obj.mode = mode;
            obj.s = s;
        end

        function obj = fuzzyfication(obj)
            if strcmp(obj.s, 'linear')
                obj.R = cell(1,5);
                sigma = [12 15 15 15 12];
                mean = [-40 -20 0 20 40];
                for i = 1:length(obj.R)
                    obj.R{i} = gaussmf(obj.u, [sigma(i) mean(i)]);
                end
            else
                obj.R = cell(1,3);
                sigma = [20 20 20];
                mean = [-25 0 25];
                for i = 1:length(obj.R)
                    obj.R{i} = gaussmf(obj.u, [sigma(i) mean(i)]);
                end
            end

            a = ones(size(obj.R));
            b = ones(size(obj.R));

            % MNK: Minimalizacja sumy kwadratów błędów
            fun = @(params) sum((obj.y_ucz - obj.fuzzy_system(params, obj.u_ucz, obj.u_start, obj.u_end, obj.R)).^2);
            initial_params = [a(1), b(1)];
            for i = 2:length(obj.R)
                initial_params = [initial_params, a(i), b(i)];
            end
            opt = optimset('MaxFunEvals', 10^6, 'MaxIter', 10^6);
            obj.params = fminsearch(fun, initial_params, opt);

            obj.y_mod_ucz = obj.fuzzy_system(obj.params, obj.u_ucz, obj.u_start, obj.u_end, obj.R);
            obj.y_mod_wer = obj.fuzzy_system(obj.params, obj.u_wer, obj.u_start, obj.u_end, obj.R);

            E_ucz = sum((obj.y_ucz - obj.y_mod_ucz).^2)/(obj.n/2);
            E_wer = sum((obj.y_wer - obj.y_mod_wer).^2)/(obj.n/2);
            fprintf('Wartość E_ucz = %.3f (F_1)\n', E_ucz);
            fprintf('Wartość E_wer = %.3f (F_1)\n', E_wer);
        end

        function obj = fuzzyfication_disturbance(obj)
            if strcmp(obj.s, 'linear')
                obj.R_disturbance = cell(1,5);
                sigma = [10 12 15 12 10];
                mean = [-15 -10 0 10 15];
                for i = 1:length(obj.R_disturbance)
                    obj.R_disturbance{i} = gaussmf(obj.u_disturbance, [sigma(i) mean(i)]);
                end
            else
                obj.R_disturbance = cell(1,3);
                sigma = [15 15 15];
                mean = [-15 0 15];
                for i = 1:length(obj.R_disturbance)
                    obj.R_disturbance{i} = gaussmf(obj.u_disturbance, [sigma(i) mean(i)]);
                end
            end

            a = ones(size(obj.R_disturbance));
            b = ones(size(obj.R_disturbance));

            % MNK: Minimalizacja sumy kwadratów błędów
            fun = @(params) sum((obj.y_ucz_disturbance - obj.fuzzy_system(params, obj.u_ucz_disturbance, obj.u_start_disturbance, obj.u_end_disturbance, obj.R_disturbance)).^2);
            initial_params = [a(1), b(1)];
            for i = 2:length(obj.R_disturbance)
                initial_params = [initial_params, a(i), b(i)];
            end
            opt = optimset('MaxFunEvals', 10^6, 'MaxIter', 10^6);
            obj.params_disturbance = fminsearch(fun, initial_params, opt);
            
            obj.y_mod_ucz_disturbance = obj.fuzzy_system(obj.params_disturbance, obj.u_ucz_disturbance, obj.u_start_disturbance, obj.u_end_disturbance, obj.R_disturbance);
            obj.y_mod_wer_disturbance = obj.fuzzy_system(obj.params_disturbance, obj.u_wer_disturbance, obj.u_start_disturbance, obj.u_end_disturbance, obj.R_disturbance);

            E_ucz = sum((obj.y_ucz_disturbance - obj.y_mod_ucz_disturbance).^2)/(obj.n/2);
            E_wer = sum((obj.y_wer_disturbance - obj.y_mod_wer_disturbance).^2)/(obj.n/2);
            fprintf('Wartość E_ucz = %.3f (F_D)\n', E_ucz);
            fprintf('Wartość E_wer = %.3f (F_D)\n', E_wer);
        end

        function Y = fuzzy_system(obj, params, u, u_start, u_end, R)
            a = zeros(size(R));
            b = zeros(size(R));
            for i = 1:length(R)
                a(i) = params(2*i-1);
                b(i) = params(2*i);
            end
        
            Y = zeros(size(u));
        
            for i = 1:length(u)
                index = round((u(i)-u_start)*obj.n / (u_end-u_start));

                if index <= 0
                    index = 1;
                elseif index > obj.n
                    index = obj.n;
                end

                weight = 0;
                if strcmp(obj.s, 'linear')
                    for j = 1:length(R)
                        Y(i) = Y(i) + R{j}(index)*(a(j) + b(j)*u(i));
                        weight = weight + R{j}(index);
                    end
                else
                    for j = 1:length(R)
                        % Y(i) = Y(i) + R{j}(index)*(a(j) + tanh(b(j)*u(i)));
                        Y(i) = Y(i) + R{j}(index)*(a(j)+sinh(b(j)*u(i)));
                        % Y(i) = Y(i) + R{j}(index)*(a(j)+b(j)*u(i)^2);
                        weight = weight + R{j}(index);
                    end
                end
                Y(i) = Y(i) / weight;
            end
        end
        
        function model(obj, a, b, K, u, y, kk, delay)
            y_mod = zeros(size(y));
            y_mod(1:delay+2) = y(1:delay+2);

            z = zeros(1, kk);
            z_disturbance = zeros(1, kk);
            z(1:delay+2) = u(1, 1:delay+2);
            z_disturbance(1:delay+2) = u(2, 1:delay+2);

            if strcmp(obj.mode, 'ARX')
                for k = delay+3:kk
                    index = round((u(1,k)-obj.u_start)*obj.n / (obj.u_end-obj.u_start));
                    index_disturbance = round((u(2,k)-obj.u_start_disturbance)*obj.n / (obj.u_end_disturbance-obj.u_start_disturbance));
                    
                    z(k) = obj.find_value(u(1,k), obj.R, obj.params, index);
                    z_disturbance(k) = obj.find_value(u(2,k), obj.R_disturbance, obj.params_disturbance, index_disturbance);
                    
                    y_mod(k) = - a*[y(k-1:-1:k-2)]' + b/K*[z(k-(delay+1):-1:k-(delay+2))]' + b/K*[z_disturbance(k-1:-1:k-2)]';
                end
            else
                for k = delay+3:kk
                    index = round((u(1,k)-obj.u_start)*obj.n / (obj.u_end-obj.u_start));
                    index_disturbance = round((u(2,k)-obj.u_start_disturbance)*obj.n / (obj.u_end_disturbance-obj.u_start_disturbance));
                    
                    z(k) = obj.find_value(u(1,k), obj.R, obj.params, index);
                    z_disturbance(k) = obj.find_value(u(2,k), obj.R_disturbance, obj.params_disturbance, index_disturbance);

                    y_mod(k) = - a*[y_mod(k-1:-1:k-2)]' + b/K*[z(k-(delay+1):-1:k-(delay+2))]' + b/K*[z_disturbance(k-1:-1:k-2)]';
                end
            end

            E = sum((y - y_mod).^2)/(kk);
            fprintf('Model %s Hammerstein \n', obj.mode);
            fprintf('E = %.3f \n', E);

            figure;
            hold on;
            stairs(0:kk-1, y_mod, 'b-', 'LineWidth', 1.2);
            stairs(0:kk-1, y, 'r-', 'LineWidth', 0.8);
            hold off;
            xlabel('k');
            ylabel('y(k)');
            title(sprintf('Model %s Hammerstein \n E = %.3f', obj.mode, E));
            legend('y_{mod}', 'y');
            grid on;
        end

        function Y = find_value(obj, u, R, params, index)
            a = zeros(size(R));
            b = zeros(size(R));
            for i = 1:length(R)
                a(i) = params(2*i-1);
                b(i) = params(2*i);
            end
        
            Y = 0;
            weight = 0;

            if index <= 0
                index = 1;
            elseif index > obj.n
                index = obj.n;
            end

            if strcmp(obj.s, 'linear')
                for i = 1:length(R)
                    Y = Y + R{i}(index)*(a(i) + b(i)*u);
                    weight = weight + R{i}(index);
                end
            else
                for i = 1:length(R)
                    % Y = Y + R{i}(index)*(a(i) + tanh(b(i)*u));
                    Y = Y + R{i}(index)*(a(i)+sinh(b(i)*u));
                    % Y = Y + R{i}(index)*(a(i)+b(i)*u^2);
                    weight = weight + R{i}(index);
                end
            end
            Y = Y / weight;
        end

        function show_fuzzy_system(~, u, u_ucz, u_wer, y_ucz, y_mod_ucz, y_wer, y_mod_wer, R, enforce)
            figure;
            hold on;
            for i = 1:length(R)
                plot(u, R{i});
            end
            hold off;
            grid on;
            xlabel('u');
            ylabel('y');
            xlim([u(1) u(end)]);
            title(sprintf('Zbiory rozmyte y(u) [ %s ]', enforce));

            figure;
            hold on;
            plot(u_ucz, y_ucz, 'b-');
            plot(u_ucz, y_mod_ucz, 'ro');
            hold off;
            xlim([u_ucz(1) u_ucz(end)]);
            ylim([y_ucz(1) y_ucz(end)]);
            grid on;
            xlabel('u');
            ylabel('y');
            title(sprintf('Charakterystyka statyczna y_{ucz}(u_{ucz}) [ %s ] \n Następniki liniowe', enforce));
            legend('y(u)', 'y(u) - fuzzy', 'Location', 'northwest');

            figure;
            hold on;
            plot(u_wer, y_wer, 'b-');
            plot(u_wer, y_mod_wer, 'ro');
            hold off;
            xlim([u_wer(1) u_wer(end)]);
            ylim([y_wer(1) y_wer(end)]);
            grid on;
            xlabel('u');
            ylabel('y');
            title(sprintf('Charakterystyka statyczna y_{wer}(u_{wer}) [ %s ] \n Następniki liniowe', enforce));
            legend('y(u)', 'y(u) - fuzzy', 'Location', 'northwest');
        end

        function show_static_characteristic(~, u, y, enforce)
            figure;
            plot(u, y, 'b', 'LineWidth', 1.2);
            xlabel('u');
            ylabel('y');
            xlim([u(1) u(end)]);
            ylim([y(1) y(end)]);
            legend('y(u)', 'Location', 'northwest');
            title(sprintf("Charakterystyka statyczna w zależności od %s", enforce));
            grid on;
        end
    end
end