clear all;
close(figure(1));
close(figure(2));

%================= INIT ====================%
Tp = 10; % okres próbkowania 
k = 9000/Tp; % liczba próbek

time = linspace(-200,(k-1)*Tp-200,k);

A1 = 485; % cm^2
C2 = 0.7;
a1 = 21;
a2 = 22;

F1 = zeros(1,k);
F2 = zeros(1,k);
F3 = zeros(1,k);
F1in = zeros(1, k);

F1_L = zeros(1,k);
F2_L= zeros(1,k);
F3_L = zeros(1,k);

V1 = zeros(1, k);
V2 = zeros(1, k);
dV1 = zeros(1, k);
dV2 = zeros(1, k);

V1L = zeros(1, k);
V2L = zeros(1, k);
dV1L = zeros(1, k);
dV2L = zeros(1, k);

h2 = zeros(1, k);
h1 = zeros(1, k);

h2_L = zeros(1, k);
h1_L = zeros(1, k);

tau = 200/Tp; % s

h2_0 = 17.894;
h1_0 = a2^2/a1^2*h2_0;
V1_0 = h1_0*A1;
V2_0 = C2*h2_0^3;


h2(1:tau) = h2_0;
h1(1:tau) = h1_0;
h2_L(1:tau) = h2_0;
h1_L(1:tau) = h1_0;

FD = 18;

%================================================%
%========= Zmienne liniowy model rozmyty ========%
%================================================%
    n = 5;
    h2_0_L = zeros(1, n);
    h1_0_L = zeros(1, n);
    V1_0_L = zeros(1, n);
    V1_0_L = zeros(1, n);
    
    sqrt_h1_0 = zeros(1, n);
    sqrt_h2_0 = zeros(1, n);
    
    sqrtV20 = zeros(1, n);
    sqrt2V20 = zeros(1, n);

    w = zeros(1,n);
      
    h1_L_M = zeros(1,n);
    h2_L_M = zeros(1,n);
    V1L_M = zeros(1,n);
    V2L_M = zeros(1,n);





%=============== Przynależności wejściowe ===============%


h2_0_L(1) = 5;
h2_0_L(2) = 12;
h2_0_L(3) = 19;
h2_0_L(4) = 26;
h2_0_L(5) = 33;

h2_max = 40;
res = 1000;
x = linspace(0, h2_max, res);
%diff = 3.25;
diff = 3.25;

% w_arr{1} = trimf(x, [0  h2_0_L(1)     h2_0_L(2)+diff]);  % Przynależność do zbioru A1
% w_arr{2} = trimf(x, [h2_0_L(2)-diff   h2_0_L(2)     h2_0_L(2)+diff]);  % Przynależność do zbioru A2
% w_arr{3} = trimf(x, [h2_0_L(3)-diff   h2_0_L(3)     h2_0_L(3)+diff]); % Przynależność do zbioru A3
% w_arr{4} = trimf(x, [h2_0_L(4)-diff   h2_0_L(4)     h2_0_L(4)+diff]);  % Przynależność do zbioru A4
% w_arr{5} = trimf(x, [h2_0_L(5)-diff   h2_0_L(5)            h2_max]);  % Przynależność do zbioru A5

w_arr{1} = gaussmf(x, [diff   h2_0_L(1)]);  % Przynależność do zbioru A1
w_arr{2} = gaussmf(x, [diff   h2_0_L(2)]);  % Przynależność do zbioru A2
w_arr{3} = gaussmf(x, [diff   h2_0_L(3)]); % Przynależność do zbioru A3
w_arr{4} = gaussmf(x, [diff   h2_0_L(4)]);  % Przynależność do zbioru A4
w_arr{5} = gaussmf(x, [diff   h2_0_L(5)]);  % Przynależność do zbioru A5

% w_arr{1} = trapmf(x, [-5 2.5  7.5 10]);  % Przynależność do zbioru A1
% w_arr{2} = trapmf(x, [5 10 14 19]);  % Przynależność do zbioru A2
% w_arr{3} = trapmf(x, [12 17 21 26]); % Przynależność do zbioru A3
% w_arr{4} = trapmf(x, [19 24 28 33]);  % Przynależność do zbioru A4
% w_arr{5} = trapmf(x, [26 31 35 40]);  % Przynależność do zbioru A5

figure(3)
plot(x, w_arr{1}, x, w_arr{2}, x, w_arr{3}, x, w_arr{4}, x, w_arr{5});
title('Przynależności wejściowe');
    

for j = 1:n
    h1_0_L(j) = a2^2/a1^2*h2_0_L(j);
    V1_0_L(j) = h1_0_L(j)*A1;
    V2_0_L(j) = C2*h2_0_L(j)^3;
    
    sqrt_h1_0(j) = sqrt(h1_0_L(j));
    sqrt_h2_0(j) = sqrt(h2_0_L(j));
    
    sqrtV20(j) = nthroot(V2_0_L(j)/C2, 3);
    sqrt2V20(j) = nthroot((V2_0_L(j)/C2)^2, 3);
end

%=======================================%
V1(1:tau) = V1_0;
V2(1:tau) = V2_0;
V1L(1:tau) = V1_0;
V2L(1:tau) = V2_0;

%===============================%
%========= F1in wartość ========%
   ster = 75/2;
   ster2 = 72*1.5;
   ster3 = 75;
%===============================%

F1in(1:tau) = 75;
F1in(tau+1:k/2) = ster;
F1in(k/4+1:k) = ster2;
F1in(k/2+1:k) = ster3;
F1(1:tau) = 75;
F1_L(1:tau) = ster;
F2(1:tau) = a1*sqrt(h1_0);
F3(1:tau) = a2*sqrt(h2_0);
F2_L(1:tau) = F2(1:tau);
F3_L(1:tau) = F3(1:tau);





for i=(tau+1):k
    %================== MODEL NIELINIOWY ==================%  
    % zmodyfikowana metoda Eulera

    pomV1 = V1(i-1) + 0.5*Tp*(F1(i-1) + FD - a1*sqrt(h1(i-1)));
    if pomV1 < 0
        pomV1 = 0;
    end
    pomV2 = V2(i-1) + 0.5*Tp*(a1*sqrt(h1(i-1)) - a2*sqrt(h2(i-1)));
    if pomV2 < 0
        pomV2 = 0;
    end

    pomh1 = pomV1/A1;
    pomh2 = nthroot(pomV2/C2, 3);

    F1(i) = F1in(i - tau);
    F2(i) = a1*sqrt(pomh1);
    F3(i) = a2*sqrt(pomh2);

    dV1(i) = F1(i) + FD - F2(i); 
    dV2(i) = F2(i) - F3(i);

    V1(i) = V1(i-1) + Tp*dV1(i);

    if V1(i) < 0
        V1(i) = 0;
    end

    V2(i) = V2(i-1) + Tp*dV2(i);

    if V2(i) < 0
        V2(i) = 0;
    end

    h1(i) = V1(i)/A1;
    h2(i) = nthroot(V2(i)/C2, 3);



  %================== LINIOWY MODEL ROZMYTY ==================% 

    index = round(res   *  min(h2_max , h2_L(i-1))    /h2_max);
    
    for j = 1:n
        w(j) = w_arr{j}(index);
    end

  for j = 1:n
        if w(j) == 0
            V1L_M(j) = 0;
            V2L_M(j) = 0;
            h1_L_M(j) = 0;
            h2_L_M(j) = 0;
            continue
        end


       sqrt_h1_L = sqrt_h1_0(j) + (h1_L(i-1) - h1_0_L(j))/(2*sqrt_h1_0(j));
       sqrt_h2_L = sqrt_h2_0(j) + (h2_L(i-1) - h2_0_L(j))/(2*sqrt_h2_0(j));
    
       pomV1L = V1L(i-1) + 0.5*Tp*(F1_L(i-1) + FD - a1*sqrt_h1_L);
        if pomV1L < 0
            pomV1L = 0;
        end
    
       pomV2L = V2L(i-1) + 0.5*Tp*(a1*sqrt_h1_L - a2*sqrt_h2_L);
       if pomV2L < 0 
           pomV2L = 0;
       end
    
       pomh1_L = pomV1L/A1;
       pomh2_L = sqrtV20(j) + (pomV2L - V2_0_L(j))/(3*sqrt2V20(j));
       sqrt_h1_L = sqrt_h1_0(j) + (pomh1_L - h1_0_L(j))/(2*sqrt_h1_0(j));
       sqrt_h2_L = sqrt_h2_0(j) + (pomh2_L - h2_0_L(j))/(2*sqrt_h2_0(j));
    
       F1_L(i) = F1in(i - tau);
       F2_L(i) = a1*sqrt_h1_L ;
       F3_L(i) = a2*sqrt_h2_L ;
    
       dV1L(i) = F1_L(i) + FD - F2_L(i); 
       dV2L(i) = F2_L(i) - F3_L(i);
    
       V1L_M(j) = V1L(i-1) + Tp*dV1L(i);
    
       if V1L_M(j) < 0
           V1L_M(j) = 0;
       end
    
       V2L_M(j) = V2L(i-1) + Tp*dV2L(i);
    
       if V2L_M(j) < 0
           V2L_M(j) = 0;
       end
    
       h1_L_M(j) = V1L_M(j)/A1;
       h2_L_M(j) = sqrtV20(j) + (V2L_M(j) - V2_0_L(j))/(3*sqrt2V20(j));
  end

%================================%

   h1_L_M_suma = 0;
   h2_L_M_suma = 0;
   V1L_M_suma = 0;
   V2L_M_suma = 0;

   for j = 1:n
        h1_L_M_suma = h1_L_M_suma + h1_L_M(j)*w(j);
        h2_L_M_suma = h2_L_M_suma + h2_L_M(j)*w(j);
        V1L_M_suma = V1L_M_suma + V1L_M(j)*w(j);
        V2L_M_suma = V2L_M_suma + V2L_M(j)*w(j);
        
   end
        w_suma = sum(w(1:n));

        h1_L(i) = h1_L_M_suma/w_suma;
        h2_L(i) = h2_L_M_suma/w_suma;
        V1L(i) = V1L_M_suma/w_suma;
        V2L(i) = V2L_M_suma/w_suma;

end



%================ WYKRESY ================%

figure(2)
hold on
plot(time, V1,'--')
title(strcat("F1in = ", int2str(ster)," cm^{3}/h"))
xlabel('czas [s]') 
ylabel('objętość [cm]') 
plot(time, V2,'--')
plot(time, V1L,Color=[0,0.75,1])
plot(time, V2L,Color=[0.7,0.2,0.2])
xlim([0 k*Tp-200])
legend("V1", "V2", "V1L", "V2L",'Location','northwest')

figure(1)
hold on
plot(time, h1,'--')
title(strcat("F1in = ", int2str(ster)," cm^{3}/h"))
xlabel('czas [s]') 
ylabel('wysokość [cm]') 
%hold on
plot(time, h2,'--')

plot(time, h1_L,Color=[0,0.75,1])
plot(time, h2_L,Color=[0.7,0.2,0.2])
xlim([0 k*Tp-200])
legend("h1", "h2", "h1L", "h2L",'Location','northwest')
hold off




hold off


