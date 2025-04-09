%%
% Wczytanie danych
data = readmatrix('samplingData_11_02_2025_11_59_44.txt');
Y1 = data(:,2); % Druga kolumna
Y2 = data(:,3);
t = 1:5:5*length(Y1);
plot(t, Y1);
hold on;
plot(t, Y2);
legend('Temperatura komory', 'Temperatura grzałki', 'Location', 'best');
xlabel('t[s]');
ylabel('[C]');
grid on;

saveas(gcf, 'D:/EiTI/MGR/raporty/raport_luty/pictures/step_response.png');  % Zapisuje jako plik PNG

%%


% Pobranie drugiej i trzeciej kolumny
x = (1:length(data(:,2))-500)';  % Indeksy czasowe
y1 = data(1:end-500,2); % Druga kolumna
y2 = data(1:end-500,3); % Trzecia kolumna

% --- Regresja wielomianowa i ekstrapolacja ---
poly_order = 5; % Stopień wielomianu (można dostosować)
p1 = polyfit(x, y1, poly_order); % Dopasowanie do drugiej kolumny
p2 = polyfit(x, y2, poly_order); % Dopasowanie do trzeciej kolumny

% Tworzymy nowe punkty x do ekstrapolacji
x_new = (1:max(x) + 500)'; % Wydłużamy o 200 punktów
y1_extrap = polyval(p1, x_new);
y2_extrap = polyval(p2, x);

% --- Filtr Savitzky-Golay dla wygładzenia ---
window_size = 11; % Długość okna
poly_order_sg = 3; % Stopień wielomianu dla filtru

y1_smooth = sgolayfilt(y1_extrap, poly_order_sg, window_size);
y2_smooth = sgolayfilt(y2_extrap, poly_order_sg, window_size);

% --- Znalezienie miejsca wypłaszczenia ---
tolerance = 0.001; % Próg zmian wartości do uznania za wypłaszczenie

idx1 = find(abs(diff(y1_smooth)) < tolerance, 1); % Pierwszy punkt stabilizacji
idx2 = find(abs(diff(y2_smooth)) < tolerance, 1);

flat_point_x1 = x_new(idx1);
flat_point_x2 = x_new(idx2);

flat_point_y1 = y1_smooth(idx1);
flat_point_y2 = y2_smooth(idx2);

% --- Wizualizacja wyników ---
figure;
subplot(2,1,1);
plot(x, y1, 'r--', 'LineWidth', 1.2); hold on;
plot(x_new, y1_extrap, 'g-', 'LineWidth', 1.5); % Regresja wielomianowa + ekstrapolacja
plot(x_new, y1_smooth, 'b-', 'LineWidth', 1.5); % Savitzky-Golay
plot(flat_point_x1, flat_point_y1, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k'); % Punkt wypłaszczenia
legend('Oryginalne dane', 'Ekstrapolacja', 'Savitzky-Golay', 'Punkt wypłaszczenia');
title('Druga kolumna - Aproksymacja i Wypłaszczenie');
xlabel('Indeks czasu');
ylabel('Wartość');
grid on;

subplot(2,1,2);
plot(x, y2, 'r--', 'LineWidth', 1.2); hold on;
plot(x, y2_extrap, 'g-', 'LineWidth', 1.5);
plot(x, y2_smooth, 'b-', 'LineWidth', 1.5);
plot(flat_point_x2, flat_point_y2, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
legend('Oryginalne dane', 'Ekstrapolacja', 'Savitzky-Golay', 'Punkt wypłaszczenia');
title('Trzecia kolumna - Aproksymacja i Wypłaszczenie');
xlabel('Indeks czasu');
ylabel('Wartość');
grid on;

%%
                % Wartości początkowe dla solvera
                % x0 = [pH_1, pH_2, pH]; % Początkowe przybliżenia pH
                % x0 = pH;

                % Rozwiązywanie układu równań
                % options = optimoptions('fsolve', ...
                %         'Display', 'off', ... % Wyświetla kolejne kroki iteracji
                %         'TolFun', 1e-6, ...   % Zmniejsza tolerancję funkcji (mniejsze błędy)
                %         'TolX', 1e-12, ...     % Zmniejsza tolerancję zmiennej (dokładniejsze `x`)
                %         'MaxIterations', 1000, ... % Zwiększa liczbę iteracji
                %         'Algorithm', 'levenberg-marquardt'); % Alternatywny algorytm
                % [sol, fval, exitflag]= fsolve(@(x) obj.calculate_pH(x, C_1e, C_2e), x0, options);
                % if (exitflag == -1)
                %     disp(['Exit flag: ', num2str(exitflag)]);
                %     disp(['Residuum: ', num2str(norm(fval))]);
                % end
                % 
                % % Zapisanie wyników
                % pH_1 = sol(1);
                % pH_2 = sol(2);
                % pH = sol(3);
                % [pH_1, pH_2, pH] = obj.calculate_pH(C_1e, C_2e, C_1in(i), C_2in(i), pH_1, pH_2, pH);

                % C_1in_curr = C_1in(i);
                % C_2in_curr = C_2in(i);
                % C_1e_curr = C_1e;
                % C_2e_curr = C_2e;
                
                % % Symboliczne równania
                % if (i == 2)
                %     % eq1 = 10^(-3*x1) + 10^(-2*x1) * obj.Ka + 10^(-x1) * (-obj.Ka * C_1in_curr - obj.Kw) - obj.Ka * obj.Kw == 0;
                %     % eq2 = 10^(-3*x2) + 10^(-2*x2) * (obj.Ka + C_2in_curr) + 10^(-x2) * (obj.Ka * C_2in_curr - obj.Kw) - obj.Ka * obj.Kw == 0;
                %     eq3 = 10^(-3*x3) + 10^(-2*x3) * (obj.Ka + C_2e_curr) + 10^(-x3) * (obj.Ka * (C_2e_curr - C_1e_curr) - obj.Kw) - obj.Ka * obj.Kw == 0;
                % 
                %     % Rozwiązanie symboliczne
                %     % sol = vpasolve([eq1, eq2, eq3], [x1, x2, x3]);
                %     sol = vpasolve([eq3], [x3]);
                % 
                %     % Zamiana na wartości numeryczne
                %     % pH_1 = double(sol.x1);
                %     % pH_2 = double(sol.x2);
                %     pH = double(sol);
                % 
                %     % phFunc = matlabFunction([sol.x1, sol.x2, sol.x3], 'Vars', [C1e, C2e, C1in, C2in]);
                %     phFunc = matlabFunction([eq3], 'Vars', [C1e, C2e]);
                % else
                %     pH = phFunc(C_1e_curr, C_2e_curr);
                % end

                % eq1 = 10^(-3*x1) + 10^(-2*x1) * obj.Ka + 10^(-x1) * (-obj.Ka * C_1in(i) - obj.Kw) - obj.Ka * obj.Kw == 0;
                % eq2 = 10^(-3*x2) + 10^(-2*x2) * (obj.Ka + C_2in(i)) + 10^(-x2) * (obj.Ka * C_2in(i) - obj.Kw) - obj.Ka * obj.Kw == 0;
                % eq3 = 10^(-3*x3) + 10^(-2*x3) * (obj.Ka + C_2e) + 10^(-x3) * (obj.Ka * (C_2e - C_1e) - obj.Kw) - obj.Ka * obj.Kw == 0;

                % Rozwiązanie symboliczne
                % sol = vpasolve([eq1, eq2, eq3], [x1, x2, x3]);
                % sol = vpasolve([eq3], [x3]);

                % Zamiana na wartości numeryczne
                % pH_1 = double(sol.x1);
                % pH_2 = double(sol.x2);
                % pH = double(sol);