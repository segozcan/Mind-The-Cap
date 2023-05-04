
clear all
% Specify project requirements

Vimax = 18;
Vimin = 12;
Vo = 48;
Po = 48;

% Choose an estimated efficiency
eff_est = 0.88;
% Choose transformer current ripple ratio
ripple_ratio = 0.40;

%% Duty Range vs Turns Ratio
% We observe the effect of turns ratio on the duty range of our converter
% given the input and output voltage range

% We sweep N2/N1 between 1 and 20
N = linspace(1,20,20);

Dmax_vec = zeros(1,20);
Dmin_vec = zeros(1,20);

% Calculate min and max duty values for all N values
for i=1:20
    d1 = ((1/N(i)) * (Vo/Vimin))./(1+(1/N(i)) * (Vo/Vimin));
    d2 = ((1/N(i)) * (Vo/Vimax))./(1+(1/N(i)) * (Vo/Vimax));
    dmax = max(d1,d2);
    dmin = min(d1,d2);
    Dmax_vec(i) = dmax;
    Dmin_vec(i) = dmin;
end

figure;
plot(N,Dmax_vec);
hold on
plot(N,Dmin_vec);
hold on
legend
hold off

%% Determine Lm
clear dmin dmax Dmin_vec Dmax_vec d1 d2 i N
Vin_vec = linspace(12, 18, 50);

% Choose max duty as 0.5, calculate N
Dmax = 0.5;
N = (Vo/Vimin)*(1-Dmax)/Dmax;
% Duty range is between 0.4 & 0.5
D_vec = (Vo./(Vin_vec*N))./(1+Vo./(Vin_vec*N));

% Choose switching frequency
fs = 2e5; % 200kHz

% Determine the transformer current ripple
Pin = Po/eff_est;
I_pri_avg = Pin./Vin_vec;
Delta_I_Lm = I_pri_avg.*ripple_ratio;

Lmvec = zeros(1,50);

% Find required Lm for different input voltages
for i=1:50
    Lmvec(i) = (Vin_vec(i)*D_vec(i))/(fs*Delta_I_Lm(i));
end

figure;
plot(Vin_vec, Lmvec);

% Choose Lm as 30uH to satisfy ripple ratio requirement for all inputs
Lm = 30e-6;

%% Choose Gap Length and Turn Number

% Permeability of air
mu0 = 4*pi*1e-7;

% Define core parameters
A = 30.80e-3;
B = 19.50e-3;
C = 7.20e-3;
F = 7.30e-3;

% Calculate areas from core parameters
Asmall = F*(A-B)/2;
Alarge = C*F;

% Sweep primary turn numbers and calculate other parameters
turns_pri_vec = linspace(1, 20, 20);
reluctance_vec = (turns_pri_vec.^2)./Lm;
gap_interval_vec = reluctance_vec.*(2*mu0*Asmall*Alarge)./(Alarge+2*Asmall);

fig = figure;
plot(turns_pri_vec, gap_interval_vec)
title("Gap length vs turn number")
figure;
plot(turns_pri_vec, reluctance_vec)
grid on

%exportgraphics(fig, "deneme.pdf")

% Choose a gap value (A4 width = 0.1mm), corresponds to 9 turns
gap = 0.1e-3;
N_pri = 9;

%% Determine Cable Length and Type

