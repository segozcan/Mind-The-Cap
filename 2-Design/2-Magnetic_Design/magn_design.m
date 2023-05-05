close all
clear 
%%
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
N = (Vo/Vimin)*(1-Dmax)/Dmax; % N2/N1 = 4
% Duty range is between 0.4 & 0.5
D_vec = (Vo./(Vin_vec*N))./(1+Vo./(Vin_vec*N));

% Choose switching frequency
fs = 2e5; % 200kHz

% Determine the transformer current ripple
Pin = Po/eff_est;
I_in_avg = Pin./Vin_vec;
I_Lm_avg = I_in_avg./D_vec;
Delta_I_Lm_vec = I_Lm_avg.*ripple_ratio;

Lmvec = zeros(1,50);

% Find required Lm for different input voltages
for i=1:50
    Lmvec(i) = (Vin_vec(i)*D_vec(i))/(fs*Delta_I_Lm_vec(i));
end

figure;
plot(Vin_vec, Lmvec);

% Choose Lm as 12uH to satisfy ripple ratio requirement for all inputs
Lm = 12e-6;

%% Choose Gap Length and Turn Number
% https://www.ferroxcube.com/upload/media/product/file/Pr_ds/E30_15_7.pdf

% TO DO: Add core reluctance to the gap reluctance for further calculations

% Permeability of air
mu0 = 4*pi*1e-7;

% Define core parameters
A = 30.80e-3;
B = 19.50e-3;
C = 7.20e-3;
D = 15e-3;
F = 7.30e-3;
E = 9.3e-3;

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
reluctance = 6750000 ;
gap = 0.2723e-3;
N_pri = 9;

%% AWG
cu_resistivity = 1.68e-8;
skin_depth = sqrt(cu_resistivity/(pi*fs*mu0*1)) ; % in meters

strand_radius = pi*skin_depth^2*1e6 ; % 0.0668 mm^2 nearly 29 AWG

risk_factor = 0.75;
strand_pri = max(I_in_avg)/(0.182*risk_factor); % 29 AWG current rating 0.182 A
strand_sec = (Po/Vo)/(0.182*risk_factor) ; % 29 AWG current rating 0.182 A
% Determine Cable Length and Type

% Set cable diameters:

primary_parallel = strand_pri; % how many parallel cables there are 
primary_cable_diameter =  0.28702e-3;
primary_cable_count = N_pri*primary_parallel;

secondary_parallel = strand_sec; % how many parallel cables there are 
secondary_cable_diameter = 0.28702e-3;
secondary_cable_count = N_pri*N*secondary_parallel;

total_cable_area = pi*(primary_cable_diameter/2)^2*primary_cable_count +  pi*(secondary_cable_diameter/2)^2*secondary_cable_count;

window_area = ((B-C)/2)*2*E;
fill_factor_awg = total_cable_area/window_area;



%% Determine Cable Length and Type
% Set cable diameters:

primary_parallel = 2; % how many parallel cables there are 
primary_cable_diameter = 0.9e-3;
primary_cable_count = N_pri*primary_parallel;

secondary_parallel = 1; % how many parallel cables there are 
secondary_cable_diameter = 0.9e-3;
secondary_cable_count = N_pri*N*secondary_parallel;

% https://www.elektrisola.com/en/Litz-Wire/Info
litz_packing_factor = 1.28 ;
total_cable_area = pi*(primary_cable_diameter/2)^2*primary_cable_count +  pi*(secondary_cable_diameter/2)^2*secondary_cable_count;
total_cable_area = litz_packing_factor*total_cable_area;

window_area = ((B-C)/2)*2*E;
fill_factor_litz = total_cable_area/window_area;



%% Determine Max B
I_pri_peak = I_Lm_avg+I_Lm_avg.*ripple_ratio./2;
Phi = (N_pri.*(I_pri_peak))./reluctance;
B_field_sides = (Phi)./((A-B)*F); % side gap
B_field_center = (Phi)./(C*F); % center gap


plot(D_vec, B_field_center)
title("Core flux density vs duty cycle")
hold on
plot(D_vec,B_field_sides)
grid minor

%% Power Loss

% Core Loss
% https://elnamagnetics.com/wp-content/uploads/library/Ferroxcube-Materials/3C94_Material_Specification.pdf

% Coefficients (with 95% confidence bounds):
p1 =   8.892e-05;
p2 =     0.01501;
p3 =      -0.775;
p4 =       17.63;

P_density_center = p1.*(B_field_center.*1000.*ripple_ratio).^3 + p2.*(B_field_center.*1000.*ripple_ratio).^2 + p3.*(B_field_center.*1000.*ripple_ratio) + p4; % mT -> kW per m^3
P_center = P_density_center.*(C*D*F*2)*1000 ;% core center volume

P_density_sides = p1.*(B_field_sides.*1000.*ripple_ratio).^3 + p2.*(B_field_sides.*1000.*ripple_ratio).^2 + p3.*(B_field_sides.*1000.*ripple_ratio) + p4; % mT -> kW per m^3
P_sides = P_density_sides.*(((A-B)*F*2*D)+((B-C)*(D-E)*F*2))*1000 ;  % core center volume

P_core = (P_center + P_sides)*2 ; % an error margin of 200 % is given

% Copper Loss for litz wire
% AC losses proximity and skin effect can be ignored
% Length of the wires is estimated from the middle radius

lenght_pri = primary_cable_count*(2*C*1.28+2*F*1.28);
lenght_sec = secondary_cable_count*(((A+B)/2*2)+(5e-3+F)*2) ;

R_pri = (cu_resistivity*lenght_pri)/(primary_parallel*pi*(primary_cable_diameter/2)^2*primary_parallel/litz_packing_factor);
R_sec = (cu_resistivity*lenght_sec)/(secondary_parallel*pi*(secondary_cable_diameter/2)^2*secondary_parallel/litz_packing_factor);

P_cu_pri = R_pri*5.3*5.3;
P_cu_sec = R_sec*1.45*1.45;

P_cu_pri + P_cu_sec ;



