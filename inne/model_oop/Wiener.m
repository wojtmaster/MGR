classdef Wiener < Base
    properties
        v
        v_ucz
        v_wer
        v_start
        v_end
        y
        y_ucz
        y_mod_ucz
        y_wer
        y_mod_wer
        R
        params

        v_disturbance
        v_ucz_disturbance
        v_wer_disturbance
        v_start_disturbance
        v_end_disturbance
        y_disturbance
        y_ucz_disturbance
        y_mod_ucz_disturbance
        y_wer_disturbance
        y_mod_wer_disturbance
        R_disturbance
        params_disturbance

        n
        mode 
        s
    end

    methods 
        function obj = Wiener(a, b, K, n, delay, rk4, mode, s)
            u = [linspace(-60, 60, n); zeros(1, n)];
            [y, ~] = rk4(u, n);

            obj.v(1:delay+2) = y(1:delay+2);

            if strcmp(obj.mode, 'ARX')
                for k = delay+3:n
                    obj.v(k) = - a*[y(k-1:-1:k-2)]' + b/K*[u(1,k-(delay+1):-1:k-(delay+2))]';
                end
            else
                for k = delay+3:n
                    obj.v(k) = - a*[obj.v(k-1:-1:k-2)]' + b/K*[u(1,k-(delay+1):-1:k-(delay+2))]';
                end
            end

            obj.y = y(n/5:end);
            obj.y_ucz = obj.y(1:2:end-1);
            obj.y_wer = obj.y(2:2:end);
            obj.v = obj.v(n/5:end);
            obj.v_ucz = obj.v(1:2:end-1);
            obj.v_wer = obj.v(2:2:end);
            obj.v_start = obj.v(1);
            obj.v_end = obj.v(end);

            u = [zeros(1, n); linspace(-20, 20, n)];
            [y, ~] = rk4(u, n);

            obj.v_disturbance(1:delay+2) = y(1:delay+2);

            if strcmp(obj.mode, 'ARX')
                for k = delay+3:n
                    obj.v_disturbance(k) = - a*[y(k-1:-1:k-2)]' + b/K*[u(2, k-1:-1:k-2)]';
                end
            else
                for k = delay+3:n
                    obj.v_disturbance(k) = - a*[obj.v_disturbance(k-1:-1:k-2)]' + b/K*[u(2, k-1:-1:k-2)]';
                end
            end
            
            obj.y_disturbance = y(n/5:end);
            obj.y_ucz_disturbance = obj.y_disturbance(1:2:end-1);
            obj.y_wer_disturbance = obj.y_disturbance(2:2:end);
            obj.v_disturbance = obj.v_disturbance(n/5:end);
            obj.v_ucz_disturbance = obj.v_disturbance(1:2:end-1);
            obj.v_wer_disturbance = obj.v_disturbance(2:2:end);
            obj.v_start_disturbance = obj.v_disturbance(1);
            obj.v_end_disturbance = obj.v_disturbance(end);

            obj.n = 0.8*n;
            obj.mode = mode;
            obj.s = s;
        end

        function [obj] = fuzzyfication(obj)

            if strcmp(obj.s, 'linear')
                obj.R = cell(1,5);
                sigma = [12 15 15 15 12];
                mean = [-40 -20 0 20 40];
                for i = 1:length(obj.R)
                    obj.R{i} = gaussmf(obj.v, [sigma(i) mean(i)]);
                end
            else
                obj.R = cell(1,3);
                sigma = [20 20 20];
                mean = [-25 0 25];
                for i = 1:length(obj.R)
                    obj.R{i} = gaussmf(obj.v, [sigma(i) mean(i)]);
                end
            end

            a = ones(size(obj.R));
            b = ones(size(obj.R));

            fun = @(params) sum((obj.y_ucz - obj.fuzzy_system(params, obj.v_ucz, obj.v_start, obj.v_end, obj.R)).^2);
            initial_params = [a(1), b(1)];
            for i = 2:length(obj.R)
                initial_params = [initial_params, a(i), b(i)];
            end
            opt = optimset('MaxFunEvals', 10^6, 'MaxIter', 10^6);
            obj.params = fminsearch(fun, initial_params, opt);

            obj.y_mod_ucz = obj.fuzzy_system(obj.params, obj.v_ucz, obj.v_start, obj.v_end, obj.R);
            obj.y_mod_wer = obj.fuzzy_system(obj.params, obj.v_wer, obj.v_start, obj.v_end, obj.R);

            E_ucz = sum((obj.y_ucz - obj.y_mod_ucz).^2)/(obj.n/2);
            E_wer = sum((obj.y_wer - obj.y_mod_wer).^2)/(obj.n/2);
            fprintf('Wartość E_ucz = %.3f (F_1)\n', E_ucz);
            fprintf('Wartość E_wer = %.3f (F_1)\n', E_wer);
        end

        function [obj] = fuzzyfication_disturbance(obj)

            if strcmp(obj.s, 'linear')
                obj.R_disturbance = cell(1,5);
                sigma = [10 10 10 10 10];
                mean = [-15 -10 0 10 15];
                for i = 1:length(obj.R_disturbance)
                    obj.R_disturbance{i} = gaussmf(obj.v_disturbance, [sigma(i) mean(i)]);
                end
            else
                obj.R_disturbance = cell(1,3);
                sigma = [15 15 15];
                mean = [-20 0 20];
                for i = 1:length(obj.R_disturbance)
                    obj.R_disturbance{i} = gaussmf(obj.v_disturbance, [sigma(i) mean(i)]);
                end
            end

            a = ones(size(obj.R_disturbance));
            b = ones(size(obj.R_disturbance));

            fun = @(params) sum((obj.y_ucz_disturbance - obj.fuzzy_system(params, obj.v_ucz_disturbance, obj.v_start_disturbance, obj.v_end_disturbance, obj.R_disturbance)).^2);
            initial_params = [a(1), b(1)];
            for i = 2:length(obj.R_disturbance)
                initial_params = [initial_params, a(i), b(i)];
            end
            opt = optimset('MaxFunEvals', 10^6, 'MaxIter', 10^6);
            obj.params_disturbance = fminsearch(fun, initial_params, opt);

            obj.y_mod_ucz_disturbance = obj.fuzzy_system(obj.params_disturbance, obj.v_ucz_disturbance, obj.v_start_disturbance, obj.v_end_disturbance, obj.R_disturbance);
            obj.y_mod_wer_disturbance = obj.fuzzy_system(obj.params_disturbance, obj.v_wer_disturbance, obj.v_start_disturbance, obj.v_end_disturbance, obj.R_disturbance);

            E_ucz = sum((obj.y_ucz_disturbance - obj.y_mod_ucz_disturbance).^2)/(obj.n/2);
            E_wer = sum((obj.y_wer_disturbance - obj.y_mod_wer_disturbance).^2)/(obj.n/2);
            fprintf('Wartość E_ucz = %.3f (F_D)\n', E_ucz);
            fprintf('Wartość E_wer = %.3f (F_D)\n', E_wer);
        end

        function Y = fuzzy_system(obj, params, v, v_start, v_end, R)
            a = zeros(size(R));
            b = zeros(size(R));
            for i = 1:length(R)
                a(i) = params(2*i-1);
                b(i) = params(2*i);
            end
        
            Y = zeros(size(v));
        
            for i = 1:length(v)
                index = round((v(i)-v_start)*obj.n / (v_end-v_start));

                if index <= 0
                    index = 1;
                elseif index > obj.n
                    index = obj.n;
                end

                weight = 0;
                if strcmp(obj.s, 'linear')
                    for j = 1:length(R)
                        Y(i) = Y(i) + R{j}(index)*(a(j) + b(j)*v(i));
                        weight = weight + R{j}(index);
                    end
                else
                    for j = 1:length(R)
                        % Y(i) = Y(i) + R{j}(index)*(a(j) + tanh(b(j)*v(i)));
                        Y(i) = Y(i) + R{j}(index)*(a(j)+sinh(b(j)*v(i)));
                        % Y(i) = Y(i) + R{j}(index)*(a(j)+b(j)*v(i)^2);
                        weight = weight + R{j}(index);
                    end
                end
                Y(i) = Y(i) / weight;
            end
        end

        function model(obj, a, b, K, u, y, kk, delay)
            y_mod = zeros(size(y));
            y_mod(1:delay+2) = y(1:delay+2);

            v_dynamic = zeros(1, kk);
            v_dynamic_disturbance = zeros(1, kk);
            v_dynamic(1:delay+2) = y(1:delay+2);
            v_dynamic_disturbance(1:delay+2) = y(1:delay+2);

            if strcmp(obj.mode, 'ARX')
                for k = delay+3:kk
                    v_dynamic(k) = - a*[y(k-1:-1:k-2)]' + b/K*[u(1,k-(delay+1):-1:k-(delay+2))]';
                    v_dynamic_disturbance(k) = - a*[y(k-1:-1:k-2)]' + b/K*[u(2, k-1:-1:k-2)]';

                    index = round((v_dynamic(k)-obj.v_start)*obj.n / (obj.v_end-obj.v_start));
                    index_disturbance = round((v_dynamic_disturbance(k)-obj.v_start_disturbance)*obj.n / (obj.v_end_disturbance-obj.v_start_disturbance));
                    
                    y_mod(k) = obj.find_value(v_dynamic(k), obj.R, obj.params, index) + obj.find_value(v_dynamic_disturbance(k), obj.R_disturbance, obj.params_disturbance, index_disturbance);
                end
            else
                for k = delay+3:kk
                    v_dynamic(k) = - a*[v_dynamic(k-1:-1:k-2)]' + b/K*[u(1,k-(delay+1):-1:k-(delay+2))]';
                    v_dynamic_disturbance(k) = - a*[v_dynamic_disturbance(k-1:-1:k-2)]' + b/K*[u(2, k-1:-1:k-2)]';

                    index = round((v_dynamic(k)-obj.v_start)*obj.n / (obj.v_end-obj.v_start));
                    index_disturbance = round((v_dynamic_disturbance(k)-obj.v_start_disturbance)*obj.n / (obj.v_end_disturbance-obj.v_start_disturbance));

                    y_mod(k) = obj.find_value(v_dynamic(k), obj.R, obj.params, index) + obj.find_value(v_dynamic_disturbance(k), obj.R_disturbance, obj.params_disturbance, index_disturbance);
                end
            end

            E = sum((y - y_mod).^2)/(kk);
            fprintf('Model %s Wienera \n', obj.mode);
            fprintf('E = %.3f \n', E);

            figure;
            hold on;
            stairs(0:kk-1, y_mod, 'b-', 'LineWidth', 1.2);
            stairs(0:kk-1, y, 'r-', 'LineWidth', 0.8);
            hold off;
            xlabel('k');
            ylabel('y(k)');
            title(sprintf('Model %s Wienera \n E = %.3f', obj.mode, E));
            legend('y_{mod}', 'y');
            grid on;
        end

        function Y = find_value(obj, v, R, params, index)
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
                    Y = Y + R{i}(index)*(a(i) + b(i)*v);
                    weight = weight + R{i}(index);
                end
            else
                for i = 1:length(R)
                    % Y = Y + obj.R{i}(index)*(a(i) + tanh(b(i)*u));
                    Y = Y + R{i}(index)*(a(i)+sinh(b(i)*v));
                    % Y = Y + obj.R{i}(index)*(a(i)+b(i)*u^2);
                    weight = weight + R{i}(index);
                end
            end
            Y = Y / weight;
        end

        function show_fuzzy_system(~, v, v_ucz, v_wer, y_ucz, y_mod_ucz, y_wer, y_mod_wer, R, enforce)
            figure;
            hold on;
            for i = 1:length(R)
                plot(v, R{i});
            end
            hold off;
            grid on;
            xlabel('v');
            ylabel('y');
            xlim([v(1) v(end)]);
            title(sprintf('Zbiory rozmyte y(v) [ %s ]', enforce));

            figure;
            hold on;
            plot(v_ucz, y_ucz, 'b-');
            plot(v_ucz, y_mod_ucz, 'ro');
            hold off;
            xlim([v_ucz(1) v_ucz(end)]);
            ylim([y_ucz(1) y_ucz(end)]);
            grid on;
            xlabel('v');
            ylabel('y');
            title(sprintf('Charakterystyka statyczna y_{ucz}(v_{ucz}) [ %s ] \n Następniki liniowe', enforce));
            legend('y(v)', 'y(v) - fuzzy', 'Location', 'northwest');

            figure;
            hold on;
            plot(v_wer, y_wer, 'b-');
            plot(v_wer, y_mod_wer, 'ro');
            hold off;
            xlim([v_wer(1) v_wer(end)]);
            ylim([y_wer(1) y_wer(end)]);
            grid on;
            xlabel('v');
            ylabel('y');
            title(sprintf('Charakterystyka statyczna y_{wer}(v_{wer}) [ %s ] \n Następniki liniowe', enforce));
            legend('y(v)', 'y(v) - fuzzy', 'Location', 'northwest');
        end

        function show_static_characteristic(~, v, y, enforce)
            figure;
            plot(v, y, 'b', 'LineWidth', 1.2);
            xlabel('v');
            ylabel('y');
            legend('y(v)', 'Location', 'northwest');
            title(sprintf("Charakterystyka statyczna w zależności od %s", enforce));
            grid on;
        end

    end
end