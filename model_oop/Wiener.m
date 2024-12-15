classdef Wiener
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
        n
        mode 
        s
    end

    methods 
        function obj = Wiener(u, y, a, b, K, n, tau, mode, s)
            obj.y = y;
            obj.y_ucz = y(n/10:2:end-1);
            obj.y_wer = y(n/10+1:2:end);
            obj.n = 0.9*n;
            obj.mode = mode;
            obj.s = s;

            obj.v(1:tau+2) = obj.y(1:tau+2);

            if strcmp(obj.mode, 'ARX')
                for k = tau+3:n
                    obj.v(k) = - a*[y(k-1:-1:k-2)]' + b/K*[u(1,k-(tau+1):-1:k-(tau+2))]' + b*[u(2, k-1:-1:k-2)]';
                end
            else
                for k = tau+3:n
                    obj.v(k) = - a*[obj.v(k-1:-1:k-2)]' + b/K*[u(1,k-(tau+1):-1:k-(tau+2))]' + b*[u(2, k-1:-1:k-2)]';
                end
            end
            
            obj.v = obj.v(n/10:end);
            obj.v_ucz = obj.v(1:2:end-1);
            obj.v_wer = obj.v(2+1:2:end);
            obj.v_start = obj.v(1);
            obj.v_end = obj.v(end);
        end

        function [obj] = fuzzyfication(obj)

            if strcmp(obj.s, 'linear')
                obj.R = cell(1,5);
                sigma = [12 10 10 10 12];
                mean = [-30 -15 0 15 30];
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

            fun = @(params) sum((obj.y_ucz - obj.fuzzy_system(params, obj.v_ucz)).^2);
            initial_params = [a(1), b(1)];
            for i = 2:length(obj.R)
                initial_params = [initial_params, a(i), b(i)];
            end
            opt = optimset('MaxFunEvals', 10^6, 'MaxIter', 10^6);
            obj.params = fminsearch(fun, initial_params, opt);

            obj.y_mod_ucz = obj.fuzzy_system(obj.params, obj.v_ucz);
            obj.y_mod_wer = obj.fuzzy_system(obj.params, obj.v_wer);

            E_ucz = sum((obj.y_ucz - obj.y_mod_ucz).^2)/(obj.n/2);
            E_wer = sum((obj.y_wer - obj.y_mod_wer).^2)/(obj.n/2);
            fprintf('Wartość E_ucz = %.3f \n', E_ucz);
            fprintf('Wartość E_wer = %.3f \n', E_wer);

        end

        function Y = fuzzy_system(obj, params, v)
            a = zeros(size(obj.R));
            b = zeros(size(obj.R));
            for i = 1:length(obj.R)
                a(i) = params(2*i-1);
                b(i) = params(2*i);
            end
        
            Y = zeros(size(v));
        
            for i = 1:length(v)
                index = round((v(i)-obj.v_start)*obj.n / (obj.v_end-obj.v_start));

                if index <= 0
                    index = 1;
                elseif index > obj.n
                    index = obj.n;
                end

                weight = 0;
                if strcmp(obj.s, 'linear')
                    for j = 1:length(obj.R)
                        Y(i) = Y(i) + obj.R{j}(index)*(a(j) + b(j)*v(i));
                        weight = weight + obj.R{j}(index);
                    end
                else
                    for j = 1:length(obj.R)
                        % Y(i) = Y(i) + obj.R{j}(index)*(a(j) + tanh(b(j)*v(i)));
                        Y(i) = Y(i) + obj.R{j}(index)*(a(j)+sinh(b(j)*v(i)));
                        % Y(i) = Y(i) + obj.R{j}(index)*(a(j)+b(j)*v(i)^2);
                        weight = weight + obj.R{j}(index);
                    end
                end
                Y(i) = Y(i) / weight;
            end
        end

        function model(obj, a, b, u, y, K, kk, tau)
            y_mod = zeros(size(y));
            y_mod(1:tau+2) = y(1:tau+2);

            v_dynamic = zeros(1, kk);
            v_dynamic(1:tau+2) = y(1, 1:tau+2);

            if strcmp(obj.mode, 'ARX')
                for k = tau+3:kk
                    v_dynamic(k) = - a*[y(k-1:-1:k-2)]' + b/K'*[u(1,k-(tau+1):-1:k-(tau+2))]' + b*[u(2, k-1:-1:k-2)]';
                    index = round((v_dynamic(k)-obj.v_start)*obj.n / (obj.v_end-obj.v_start));
                    y_mod(k) = obj.find_value(v_dynamic(k), index);
                end
            else
                for k = tau+3:kk
                    v_dynamic(k) = - a*[v_dynamic(k-1:-1:k-2)]' + b/K'*[u(1,k-(tau+1):-1:k-(tau+2))]' + b*[u(2, k-1:-1:k-2)]';
                    index = round((v_dynamic(k)-obj.v_start)*obj.n / (obj.v_end-obj.v_start));
                    y_mod(k) = obj.find_value(v_dynamic(k), index);
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
            plot_title = sprintf('Model %s Wienera \n E = %.3f', obj.mode, E);
            title(plot_title);
            legend('y_{mod}', 'y');
            grid on;
        end

        function Y = find_value(obj, v_dynamic, index)
            a = zeros(size(obj.R));
            b = zeros(size(obj.R));
            for i = 1:length(obj.R)
                a(i) = obj.params(2*i-1);
                b(i) = obj.params(2*i);
            end
        
            Y = 0;
            weight = 0;

            if index <= 0
                index = 1;
            elseif index > obj.n
                index = obj.n;
            end

            if strcmp(obj.s, 'linear')
                for i = 1:length(obj.R)
                    Y = Y + obj.R{i}(index)*(a(i) + b(i)*v_dynamic);
                    weight = weight + obj.R{i}(index);
                end
            else
                for i = 1:length(obj.R)
                    % Y = Y + obj.R{i}(index)*(a(i) + tanh(b(i)*u));
                    Y = Y + obj.R{i}(index)*(a(i)+sinh(b(i)*v_dynamic));
                    % Y = Y + obj.R{i}(index)*(a(i)+b(i)*u^2);
                    weight = weight + obj.R{i}(index);
                end
            end
            Y = Y / weight;
        end

        function show_fuzzy_system(obj)
            figure;
            hold on;
            for i = 1:length(obj.R)
                plot(obj.v, obj.R{i});
            end
            hold off;
            grid on;
            xlabel('v');
            ylabel('y');
            xlim([-45 45]);
            title('Zbiory rozmyte y(v)');

            figure;
            hold on;
            plot(obj.v_ucz, obj.y_ucz, 'b-');
            plot(obj.v_ucz, obj.y_mod_ucz, 'ro');
            hold off;
            xlim([-25 25]);
            ylim([-40 40]);
            grid on;
            xlabel('v');
            ylabel('y');
            plot_title = sprintf('Charakterystyka statyczna y_{ucz}(v_{ucz}) \n Następniki liniowe');
            title(plot_title);
            legend('y(v)', 'y(v) - fuzzy', 'Location', 'northwest');

            figure;
            hold on;
            plot(obj.v_wer, obj.y_wer, 'b-');
            plot(obj.v_wer, obj.y_mod_wer, 'ro');
            hold off;
            xlim([-25 25]);
            ylim([-40 40]);
            grid on;
            xlabel('v');
            ylabel('y');
            plot_title = sprintf('Charakterystyka statyczna y_{wer}(v_{wer}) \n Następniki liniowe');
            title(plot_title);
            legend('y(v)', 'y(v) - fuzzy', 'Location', 'northwest');
        end

        function show_static_characteristic(obj)
            figure;
            plot(obj.v, obj.y, 'b', 'LineWidth', 1.2);
            xlabel('v');
            ylabel('y');
            legend('y(v)', 'Location', 'northwest');
            title("Charakterystyka statyczna w zależności od v");
            grid on;
        end

    end
end