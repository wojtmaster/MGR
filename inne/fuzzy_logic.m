%function [FIS, mean] = fuzzy_logic(F_10, F_D0, h_20, alpha_2)
%% Inputem jest aktualna wartość wysokości słupa cieczy w drugim zbiorniku

FIS = mamfis(...
    'NumInputs', 2, 'NumInputMFs', 5, ...
    'NumOutputs', 1, 'NumOutputMFs', 6, ...
    'AddRule', 'none');

%% Parametry funkcji przynależności
mean = struct('f_1', zeros(1, 5), 'f_d', zeros(1, 5), 'h_2', zeros(1, 6));
sigma = struct('f_1', [], 'f_d', [], 'h_2', []);
f_1 = [F_10-30, F_10-15, F_10, F_10+15, F_10+30];
f_d = [F_D0-10, F_D0-5, F_D0, F_D0+5, F_D0+10];
h_2 = [h_20-20, h_20-10, h_20, h_20+10, h_20+20, h_20+30];
% h_2 = [((f_1(1)+f_d(1))/alpha_2)^2, ((f_1(2)+f_d(2))/alpha_2)^2, ((f_1(3)+f_d(3))/alpha_2)^2, ((f_1(4)+f_d(4))/alpha_2)^2, ((f_1(5)+f_d(5))/alpha_2)^2];

sigma.f_1 = 5;   %F_1
sigma.f_d = 2;   %F_D
sigma.h_2 = 4;   %h_2
mean.f_1(:) = [f_1(1), f_1(2), f_1(3), f_1(4), f_1(5)];
mean.f_d(:) = [f_d(1), f_d(2), f_d(3), f_d(4), f_d(5)];
mean.h_2(:) = [h_2(1), h_2(2), h_2(3), h_2(4), h_2(5), h_2(6)];

%% First Input - F_1
FIS.Inputs(1).Name = 'F_1';
FIS.Inputs(1).Range = [50 130];

FIS.Inputs(1).MembershipFunctions(1).Name = 'F_10-30';
FIS.Inputs(1).MembershipFunctions(1).Type = 'gaussmf';
FIS.Inputs(1).MembershipFunctions(1).Parameters = [sigma.f_1, mean.f_1(1)];

FIS.Inputs(1).MembershipFunctions(2).Name = 'F_10-15';
FIS.Inputs(1).MembershipFunctions(2).Type = 'gaussmf';
FIS.Inputs(1).MembershipFunctions(2).Parameters = [sigma.f_1, mean.f_1(2)];

FIS.Inputs(1).MembershipFunctions(3).Name = 'F_10';
FIS.Inputs(1).MembershipFunctions(3).Type = 'gaussmf';
FIS.Inputs(1).MembershipFunctions(3).Parameters = [sigma.f_1, mean.f_1(3)];

FIS.Inputs(1).MembershipFunctions(4).Name = 'F_10+15';
FIS.Inputs(1).MembershipFunctions(4).Type = 'gaussmf';
FIS.Inputs(1).MembershipFunctions(4).Parameters = [sigma.f_1, mean.f_1(4)];

FIS.Inputs(1).MembershipFunctions(5).Name = 'F_10+30';
FIS.Inputs(1).MembershipFunctions(5).Type = 'gaussmf';
FIS.Inputs(1).MembershipFunctions(5).Parameters = [sigma.f_1, mean.f_1(5)];

% figure;
% plotmf(FIS, 'input', 1);

%% Second Input - F_D
FIS.Inputs(2).Name = 'F_D';
FIS.Inputs(2).Range = [15 45];

FIS.Inputs(2).MembershipFunctions(1).Name = 'F_D0-10';
FIS.Inputs(2).MembershipFunctions(1).Type = 'gaussmf';
FIS.Inputs(2).MembershipFunctions(1).Parameters = [sigma.f_d, mean.f_d(1)];

FIS.Inputs(2).MembershipFunctions(2).Name = 'F_D0-5';
FIS.Inputs(2).MembershipFunctions(2).Type = 'gaussmf';
FIS.Inputs(2).MembershipFunctions(2).Parameters = [sigma.f_d, mean.f_d(2)];

FIS.Inputs(2).MembershipFunctions(3).Name = 'F_D0';
FIS.Inputs(2).MembershipFunctions(3).Type = 'gaussmf';
FIS.Inputs(2).MembershipFunctions(3).Parameters = [sigma.f_d, mean.f_d(3)];

FIS.Inputs(2).MembershipFunctions(4).Name = 'F_D0+5';
FIS.Inputs(2).MembershipFunctions(4).Type = 'gaussmf';
FIS.Inputs(2).MembershipFunctions(4).Parameters = [sigma.f_d, mean.f_d(4)];

FIS.Inputs(2).MembershipFunctions(5).Name = 'F_D0+10';
FIS.Inputs(2).MembershipFunctions(5).Type = 'gaussmf';
FIS.Inputs(2).MembershipFunctions(5).Parameters = [sigma.f_d, mean.f_d(5)];

% figure;
% plotmf(FIS, 'input', 2);

%% Output - h_2
FIS.Outputs(1).Name = 'H_2';
FIS.Outputs(1).Range = [10 70];

FIS.Outputs(1).MembershipFunctions(1).Name = 'H_20-20';
FIS.Outputs(1).MembershipFunctions(1).Type = 'gaussmf';
FIS.Outputs(1).MembershipFunctions(1).Parameters = [sigma.h_2, mean.h_2(1)];

FIS.Outputs(1).MembershipFunctions(2).Name = 'H_20-10';
FIS.Outputs(1).MembershipFunctions(2).Type = 'gaussmf';
FIS.Outputs(1).MembershipFunctions(2).Parameters = [sigma.h_2, mean.h_2(2)];

FIS.Outputs(1).MembershipFunctions(3).Name = 'H_20';
FIS.Outputs(1).MembershipFunctions(3).Type = 'gaussmf';
FIS.Outputs(1).MembershipFunctions(3).Parameters = [sigma.h_2, mean.h_2(3)];

FIS.Outputs(1).MembershipFunctions(4).Name = 'H_20+10';
FIS.Outputs(1).MembershipFunctions(4).Type = 'gaussmf';
FIS.Outputs(1).MembershipFunctions(4).Parameters = [sigma.h_2, mean.h_2(4)];

FIS.Outputs(1).MembershipFunctions(5).Name = 'H_20+20';
FIS.Outputs(1).MembershipFunctions(5).Type = 'gaussmf';
FIS.Outputs(1).MembershipFunctions(5).Parameters = [sigma.h_2, mean.h_2(5)];

FIS.Outputs(1).MembershipFunctions(6).Name = 'H_20+30';
FIS.Outputs(1).MembershipFunctions(6).Type = 'gaussmf';
FIS.Outputs(1).MembershipFunctions(6).Parameters = [sigma.h_2, mean.h_2(6)];

figure;
plotmf(FIS, 'output', 1);

FIS = addRule(FIS, [ ...
    "F_1==F_10-30 & F_D==F_D0-10 => H_2=H_20-20" ...
    "F_1==F_10-30 & F_D==F_D0-5 => H_2=H_20-20" ...
    "F_1==F_10-30 & F_D==F_D0 => H_2=H_20-20" ...
    "F_1==F_10-30 & F_D==F_D0+5 => H_2=H_20-10" ...
    "F_1==F_10-30 & F_D==F_D0+10 => H_2=H_20-10" ...

    "F_1==F_10-15 & F_D==F_D0-10 => H_2=H_20-20" ...
    "F_1==F_10-15 & F_D==F_D0-5 => H_2=H_20-10" ...
    "F_1==F_10-15 & F_D==F_D0 => H_2=H_20-10" ... 
    "F_1==F_10-15 & F_D==F_D0+5 => H_2=H_20-10" ...
    "F_1==F_10-15 & F_D==F_D0+10 => H_2=H_20" ...

    "F_1==F_10 & F_D==F_D0-10 => H_2=H_20-10" ...
    "F_1==F_10 & F_D==F_D0-5 => H_2=H_20" ...
    "F_1==F_10 & F_D==F_D0 => H_2=H_20" ...
    "F_1==F_10 & F_D==F_D0+5 => H_2=H_20" ...
    "F_1==F_10 & F_D==F_D0+10 => H_2=H_20+10" ...

    "F_1==F_10+15 & F_D==F_D0-10 => H_2=H_20" ...
    "F_1==F_10+15 & F_D==F_D0-5 => H_2=H_20+10" ...
    "F_1==F_10+15 & F_D==F_D0 => H_2=H_20+10" ...
    "F_1==F_10+15 & F_D==F_D0+5 => H_2=H_20+10" ...
    "F_1==F_10+15 & F_D==F_D0+10 => H_2=H_20+20" ...

    "F_1==F_10+30 & F_D==F_D0-10 => H_2=H_20+10" ...
    "F_1==F_10+30 & F_D==F_D0-5 => H_2=H_20+20" ...
    "F_1==F_10+30 & F_D==F_D0 => H_2=H_20+20" ...
    "F_1==F_10+30 & F_D==F_D0+5 => H_2=H_20+20" ...
    "F_1==F_10+30 & F_D==F_D0+10 => H_2=H_20+30"
    ]);
%end

%% Generate input-output data and plot it
kk = 180;
alpha_2 = 20;
u = linspace(45, 135, kk)';
F_D = 30;
y = ((u+F_D)/alpha_2).^2;
% Plot of parabola
plot(u,y)
grid on
xlabel('u');
ylabel('y');
title('Nonlinear characteristics')
% Store data in appropriate form for genfis1 and anfis and plot it
data=[u y];
trndata=data(1:2:size(u),:);
chkdata=data(2:2:size(u),:);
% Plot of training and checking data generated from parabolic equation
figure;
plot(trndata(:,1), trndata(:,2), 'o', chkdata(:,1), chkdata(:,2),'x');
xlabel('u');
ylabel('y');
title('Measurement data'); 
grid on
%Initialize the fuzzy system with command genfis1. Use 5 bellshaped
%membership functions.
nu=5; 
mftype='gaussmf'; 
% fismat = genfis1(trndata, nu, mftype);
opt = genfisOptions('GridPartition');
opt.NumMembershipFunctions = nu;
opt.InputMembershipFunctionType = mftype;
fismat = genfis(trndata(:,1), trndata(:,2), opt);
%The initial membership functions produced by genfis1 are plotted
figure;
plotmf(fismat,'input',1)
xlabel('x');
ylabel('output');
title('Initial membership functions');
grid on
% Apply anfis-command to find the best FIS system - max number of
% iterations = 100;
numep = 100;
options = anfisOptions;
options.InitialFIS = fismat;
options.EpochNumber = numep;
options.ValidationData = chkdata;
% [parab, trnerr,ss,parabcheck,chkerr]=anfis(trndata,fismat,numep,[],chkdata);
[parab, trnerr, ss, parabcheck, chkerr]=anfis(trndata, options);
%Evaluate the output of FIS system using input x
anfi=evalfis(parab,u);
% Plot of trained fuzzy system using trained data
figure;
plot(trndata(:,1), trndata(:,2), 'o', chkdata(:,1), chkdata(:,2), 'x', u, anfi, '--');
grid on
xlabel('x');
ylabel('output');
title('Goodness of fit');

%save('DANE.mat', 'x', 'anfi');

%%
% Definiowanie parametrów funkcji Gaussa
sigma = 10; % Odchylenie standardowe
mean = 45; % Średnia

% Tworzenie funkcji Gaussa
x = 45:0.1:135; % Zakres wartości x
y = gaussmf(x, [sigma mean]);

% Wykres funkcji Gaussa
plot(x, y);
hold on;
title('Funkcja Gaussa');
xlabel('Wartość x');
ylabel('Przynależność');
xlim([45 135]);