function [delta_h] = set_value(kk, h_0)
    h_zad = zeros(1,kk);
    h_zad(1 : 0.25*kk) = h_0;
    h_zad(0.25*kk+1 : 0.5*kk) = h_0+15;
    h_zad(0.5*kk+1 : 0.75*kk) = h_0+20;
    h_zad(0.75*kk+1 : end) = h_0-5;
    
    delta_h = h_zad - h_0;
end