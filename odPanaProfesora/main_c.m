%% MGR
% Autor: Wojciech Rogalski
% Data: 15.12.2024r.
% Tytu�: Por�wnanie modeli Hammersteina i Winera w regulacji predykcyjnej

load ww
U = [repelem((rand(1, obiekt.kk/300) * 30 - 15), 300);
    zeros(1, obiekt.kk);
    repelem((rand(1, obiekt.kk/700) * 30 - 15), 700)];
hammerstein.testLinearModel(U, a_h, a_pH, b_h, b_pH, obiekt, 'pH', 1);