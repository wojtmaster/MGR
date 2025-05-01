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

%% Rysowanie wykresów
% u = [linspace(-45, 45, 250)
%     linspace(-15, 15, 250)*0];
% y = obiekt.static_charakteristic(u);

% u = [ones(1, obiekt.kk)*0
%     ones(1, obiekt.kk)*0];
% [y, y_L, E] = obiekt.rk4(u, obiekt.kk);
% 
% if ~exist('y_figure', 'var') || ~isvalid(y_figure)
%     y_figure = figure;
% else
%     figure(y_figure); % Jeśli istnieje, przełącz na nią
% end
% 
% plot(0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp, y, 'b-', 'LineWidth', 2);
% hold on;
% plot(0:obiekt.Tp:(obiekt.kk-1)*obiekt.Tp, y_L, 'r-', 'LineWidth', 2);
% xlabel('t [s]');
% ylabel('h_2 [cm]');
% title('Wartość sygnału wyjściowego h_2 - wymuszenie F_D');
% legend('h_2 (nieliniowe)', 'h_2 (liniowe)', 'Location', 'northeast');
% grid on;

% file = fopen('errors.txt', 'a'); % Otwórz plik do zapisu (tryb 'w' nadpisuje plik)
% fprintf(file, '%d\t%.3f\n', u(2,1), E);
% fclose(file); % Zamknij plik

% plot(u(1,:) + obiekt.F_10, y, 'b-', 'LineWidth', 2);
% title('Charakterystyka statyczna h_2 (F_1)');
% xlabel('F_1 [cm^3/s]');
% ylabel('h_2 [cm]');
% xlim([45 135]);
% grid on;

% saveas(gcf, 'D:/EiTI/MGR/raporty/raport_MGR/pictures/wymuszenie_FD.png');  % Zapisuje jako plik PNG