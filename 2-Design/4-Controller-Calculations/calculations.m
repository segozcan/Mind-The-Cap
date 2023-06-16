%% CALCULATIONS FOR UC3843A 
clc
clear all
%%
V_in = 12 ; % 18
V_out = 48 ; 
P_out = 48 ;
f_sw = 200e3 ;

L_sm= 192e-6 ;
L_pm= 12e-6 ;

D = 0.5 ; %%%% ??????

N_ps = 0.25;

I_Lpeak = 10 ;

%% Oscillator selection

syms R_t C_t

eqn = (1.72)/(R_t * C_t) == f_sw ;
s=solve(eqn,C_t)
subs(s,R_t,10e3) 

%% 
R_cs = 1 / (I_Lpeak*1.3) ; 
R_cs = 10e-3 ;

%% SLOPE COMP
R_s2 = 2.2e3
R_s1 = (1.7*R_s2*f_sw*(2*L_sm*N_ps))/(V_out*(1-D)*R_cs) - R_s2

%%
R_k = 4.7e3
R_i = (R_k*(V_out-2.5))/(2.5)

%%
R_B = 4.7e3
R_A = 10e3
R_C = 10e3
C_A = 1e-9

ESR = 100e-3
C_out = 47e-6
R_out = V_out/(P_out/V_out)

V_in_min = 12



f_p = 100e3
R_f = 10e3 
ctr = 1 ;

syms f
s_f = 2*pi*1i*f
G_opto = (R_C*ctr)/(R_f*(s_f)/(2*pi*f_p)+1)

G_BC = (R_A)/(R_B*((s_f*R_A*C_A)+1))

S_N = V_in_min*R_cs/L_pm
S_E = 1.7*(R_s2*f_sw)/(R_s1+R_s2)
Q = 1/(pi*((1+((S_E)/(S_N)))*(1-D)-0.5))

G_CO = N_ps*((1-D)/(1+D))*((s_f*ESR*C_out+1)/(s_f*R_out*C_out+1))*(1-((s_f*L_sm*D)/(R_out*(1-D)^2)))*((1/3)/(1+(s_f)/(2*pi*f_sw*Q/2)+(s_f/(2*pi*f_sw/2))^2))

F_rhpz = ((N_ps)^2)/((2*pi*L_pm*D)/(R_out*(1-D)^2))

F_c = F_rhpz/2
F_c = 10e3

G_o = 10 %%% ???????????????????????????????????????????????????

R_z = (R_i)/abs((subs(G_opto,f,F_c/5)*subs(G_BC,f,F_c/5)*G_o*subs(G_CO,f,F_c/5)))

C_z = 1/(2*pi*F_c/5*R_z)

C_p = C_z/10







