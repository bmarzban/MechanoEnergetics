% This code is the main driver for the cardiovascular mechanics model 
clear; 
tic;
flag_plot_figure = 1;
flag_swap_metabolite = 0;

%% chooosing the rat number (Mean Sham rat is Rat number 9,Mean TAC rat is number 19 )
for rat_number = 18: 19
% rat_number = 9
if rat_number<=9
    shamRat = 1;
    delta_p = 0;
else
    shamRat = 0;
end

%(Tune the following variables to fit the EDLV, ESLV, EDRV, ESLV, CO, ATP consumption Rate predicted by Ox-Phos and Mechanic model)

%         and   x_ATPase ==> ATP hydrolsis rate for LV, Septal, and RV independantly  

% adjvar = [Reference area LV , Reference area Septal,  Reference area RV,  kpassive, ksr, kforce, R_SA, R_TAC, ATP_tune_Coeff, V_LV, V_RV]

%% para set 4

adjvar_all_rest =[ 0.975 1.02 0.96 1.315 0.91 0.91 1.005 1 1.266 0.6 0.5; % rat 1 SHAM
              1.0523 0.9975 1.064 1.1 1.03 1.03 0.91 1.0 1.4265 0.6 0.5; % rat 2 SHAM
              1.095 1.05 1.02  0.14 1.63 1.63 0.72 1 1.5269 0.6 0.5; % rat 3 SHAM
              1.0486 1.0593 1.0486 1.90 0.72 0.72 1.26 1.0 1.09 0.6 0.5; % rat 4 SHAM
              1.11 1.05 0.95 1.3 0.85 0.85 0.96 1 1.436 0.6 0.5; % rat 5 SHAM
              0.995 0.99 0.96 1.47 0.91 0.91 1.03 1 1.4428 0.6 0.5; % rat 6 SHAM
              1 1.04 1.03 1.14 0.725 0.725 1.08 1 1.207 0.6 0.5; % rat 7 SHAM
              0.9 0.96 0.9 1.91 0.98 0.98 1.24 1 1.3525 0.6 0.5; % rat 8 SHAM
              1.0 1 1.0 1.09 0.83 0.83 1.0 1 1.327 0.6 0.5; % rat 9 Mean SHAM
              1.29 1.01 1.5 3.0386 0.9240 0.9240 2.2 0.46 0.695 0.6 0.6; % rat 10 TAC 1
              1.2855 0.9205 1.5152 2.1436 0.9995 0.9995 1.32 0.49 0.984 0.6 0.5; % rat 11 TAC 2
              1.3098 0.8585 0.9750 1.4568 1.3692 1.3692 1.02 0.55 1.252 0.6 0.5; % rat 12 TAC 3
              1.08 0.97 1.0 1.1676 2.5193 2.5193 1.08 0.55 1.461 0.6 0.5; % rat 13 TAC 4
              1.4128 0.94 1.46 2.2061 1.1515 1.1515 1.23 0.5 0.967 0.6 0.5; % rat 14 TAC 5
              1.105 0.93 1.34 2.3101 0.891 0.891 2.28 0.54 0.861 0.6 0.5; % rat 15 TAC 6
              1.25 0.97 1.04 1.2071 1.4 1.40 0.93 0.5 1.437 0.6 0.5; % rat 16 TAC 7
              1.357 0.95 1.56 2.1644 1.2734 1.2734 1.231 0.6 1.2167 0.6 0.5; % rat 17 TAC 8
              1.33    0.94    1.38  1.9947 1.27 1.27 1.15 0.623 1.465 0.6 0.6; % rat 18 TAC 9
              1.29 0.94 1.24 1.3494 1.34 1.34 0.95 0.6 1.529 0.6 0.5; % rat 19 TAC 10
              1.0 1 1.0 1.09 0.83 0.83 1.0 1 1.327 0.6 0.5]; % rat 20 Mean TAC
              
adjvar_all_swap = [1.0 1 1.0 1.09 0.83 0.83 1.0 1 1.327 0.6 0.5; % rat 1 SHAM - SWAP Metabolite with mean TAC
              1.0 1 1.0 1.09 0.83 0.83 1.0 1 1.327 0.6 0.5; % rat 2 SHAM - SWAP Metabolite with mean TAC
              1.0 1 1.0 1.09 0.83 0.83 1.0 1 1.327 0.6 0.5; % rat 3 SHAM - SWAP Metabolite with mean TAC
              1.0 1 1.0 1.09 0.83 0.83 1.0 1 1.327 0.6 0.5; % rat 4 SHAM - SWAP Metabolite with mean TAC
              1.0 1 1.0 1.09 0.83 0.83 1.0 1 1.327 0.6 0.5; % rat 5 SHAM - SWAP Metabolite with mean TAC
              1.0 1 1.0 1.09 0.83 0.83 1.0 1 1.327 0.6 0.5; % rat 6 SHAM - SWAP Metabolite with mean TAC
              1.0 1 1.0 1.09 0.83 0.83 1.0 1 1.327 0.6 0.5; % rat 7 SHAM - SWAP Metabolite with mean TAC
              1.0 1 1.0 1.09 0.83 0.83 1.0 1 1.327 0.6 0.5; % rat 8 SHAM - SWAP Metabolite with mean TAC
              1.0 1 1.0 1.09 0.83 0.83 1.0 1 1.327 0.6 0.5; % rat 9 Mean SHAM
              1.29 1.01 1.5 1.7 0.8250 0.8250 1 0.46 1.328 0.6 0.6; % rat 10 TAC 1 - SWAP Metabolite with mean SHAM
              1.2855 0.9205 1.5152 1.52 0.8215 0.8215 1.00 0.49 1.336 0.6 0.5; % rat 11 TAC 2 - SWAP Metabolite with mean SHAM
              1.3098 0.8585 0.975 1.25 0.8763 0.8763 1.0 0.55 1.387 0.6 0.5; % rat 12 TAC 3- SWAP Metabolite with mean SHAM
              1.08 0.97 1.0 1.18 1.0954 1.0954 1.0 0.55 1.660 0.6 0.5; % rat 13 TAC 4- SWAP Metabolite with mean SHAM
              1.4128 0.94 1.46 1.45 0.7161 0.7161 1.0 0.5 1.242 0.6 0.5; % rat 14 TAC 5- SWAP Metabolite with mean SHAM
              1.12 0.95 1.34 1.495 1 1 1 0.54 1.607 0.6 0.5; % rat 15 TAC 6- SWAP Metabolite with mean SHAM
              1.357 0.95 1.56 1.37 0.8694 0.8694 1. 0.6 1.532 0.6 0.5; % rat 16 TAC 7- SWAP Metabolite with mean SHAM
              1.25 0.97 1.04 1.3528 0.825 0.825 1 0.5 1.366 0.6 0.5; % rat 17 TAC 8- SWAP Metabolite with mean SHAM
              1.33 0.94 1.38 1.59 1.4 1.4 1. 0.623 1.719 0.6 0.6; % rat 18 TAC 9- SWAP Metabolite with mean SHAM
              1.29 0.94 1.24 1.395 0.91 0.91 1 0.6 1.5174 0.6 0.5; % rat 19 TAC 10- SWAP Metabolite with mean SHAM
              1.0 1 1.0 1.09 0.83 0.83 1.0 1 1.327 0.6 0.5]; % rat 20 Mean TAC
          
adjvar = adjvar_all_rest(rat_number,:);

if flag_swap_metabolite ==1
    adjvar = adjvar_all_swap(rat_number,:);
end


%% Read the experimental data for SHAM and TAC rats from the excel file 
data = xlsread('data1.xlsx','A3:W23');
BW  = data(rat_number , 1); % g
LVW = data(rat_number , 2); % mg
RVW = data(rat_number , 3); % mg
LW = data(rat_number , 5); % mg
HR  = data(rat_number , 6); % beats/min

edLV_target = data(rat_number , 13); % uL
esLV_target = data(rat_number , 14); % uL

SV_LV_target = edLV_target - esLV_target;
CO_target = SV_LV_target / 1000 * HR;
EF_LV_target = SV_LV_target / edLV_target * 100;
TAN = data(rat_number,16)/1000; % mole/L cell
CRtot = data(rat_number,18)/1000; % mole/L cell
% TEP = data(rat_number,20)/1000; % mole/L cell
Ox_capacity = data(rat_number,21)/data(9,21); 
Ox_capacity_sham = 1; 

% Average sham
TAN_sham = data(9,16)/1000; % mole/L cell
CRtot_sham = data(9,18)/1000; % mole/L cell

Ao = 10.260e-3; % (M per liter cell)
Co = 43.007e-3; % (M per liter cell)
Po = 35.446e-3; % (M per liter cell)
% Default settings for meant SHAM)

TEP = Po - (0.283e-3)*(Ao-TAN)/(0.082e-3); % (M per liter cell)
TEP_sham = Po - (0.283e-3)*(Ao-TAN_sham)/(0.082e-3); % (M per liter cell)
if flag_swap_metabolite ==1
    TAN = TAN_sham;
    CRtot = CRtot_sham;
    TEP = TEP_sham;
    Ox_capacity = Ox_capacity_sham;
end

if shamRat == 0
    preV = data(rat_number , 11); % mm/s
    postV = data(rat_number , 12); % mm/s
    postV = postV/1000; %m/s
    preV = preV/1000; %m/s
    rho_blood = 1060; % kg/m^3
    delta_p = 0.5*(postV^2-preV^2)*rho_blood; % Pa
    delta_p = 0.0075*delta_p; % mmHg
    if delta_p > 60
        delta_p = 31.48;
    end
    R_TAC = delta_p/CO_target*60;
else 
    R_TAC = 0;
end

CO_target = 95; % ml/min
MAP_target = 93.33; %mmHg taregt mean arterial pressure based on MAP = DBP +[1/3(SBP - DBP)];
%% Adjustable variables

R_TAC = adjvar(8)*R_TAC;

tune_ATPase_LV =  adjvar(9)* (1/ 0.6801) *1.0e-3;

Amref_LV  = adjvar(1) * 2.077 ; % LV midwall reference surface area, cm^2
Amref_SEP = adjvar(2) * Amref_LV * 0.590 ; % SEP midwall reference surface area, cm^2
Amref_RV  = adjvar(3) * 3.3 ; % RV midwall reference surface area, cm^2

% Assign initial condtion for LV and RV

% V_LV  = edLV_target/1000;% intial value for V_LV and V_RV assumed to be equal to edLV_target
% V_RV  = edLV_target/1000;%
% V_LV  = 0.6;% intial value for V_LV and V_RV assumed to be equal to edLV_target
% V_RV  = 0.4;%
V_LV  = adjvar(10);% intial value for V_LV and V_RV assumed to be equal to edLV_target
V_RV  = adjvar(11);%

Vw_LV = (LVW*2/3)/1000/1.05;
Vw_SEP =(LVW/3)/1000/1.05;
Vw_RV = RVW/1000/1.05;

%% Run the energetics model to get the metabolite concentrations

    energtics_output_LV  = EnergeticsModelScript(TAN, CRtot, TEP, Ox_capacity, tune_ATPase_LV);
   
    MgATP_LV = energtics_output_LV(1);
    MgADP_LV = energtics_output_LV(2);
    Pi_LV = energtics_output_LV(10)*1000;
    dGrATPase_LV = energtics_output_LV(6);
    
    MgATP_SEP = MgATP_LV;
    MgADP_SEP = MgADP_LV;
    Pi_SEP = Pi_LV ;
    MgATP_RV = MgATP_LV;
    MgADP_RV = MgADP_LV;
    Pi_RV = Pi_LV ;


%% Run cardiovascular mechanics model
Lsref = 1.9;
k3      = 144.5586; % transition A3 to P rate constant, 1/sec
K_T = 0.4897; 
K_D = 0.194;% Used the values from Tewari etal JMCC (9/5 BM)
alpha3  =      0.1*59.3; % Stretch sensing parameter for k3, 1/um
s3      = 9.9e-3;  % Strain at which k3 is minimum, um


%% Extract Ca coeficient based on the Ca data for diferent simulation 
para_fitted_Ca = [2	3	4	5	6	7	8	9	10;
0.0838	0.1306	0.1802	0.2557	0.3099	0.3613	0.408	0.4539	0.4602;
0.7513	0.8756	1.0274	1.4988	1.6107	1.6741	1.7902	2.1398	1.9832;
2.237	2.0486	1.948	1.852	1.6932	1.6773	1.5988	1.4952	1.4524;
0.1865	0.1815	0.1709	0.1693	0.161	0.1661	0.1425	0.1458	0.1222];
freq_all = para_fitted_Ca(1,:);
A_HR_pchip = pchip(freq_all,para_fitted_Ca(2,:));
A_HR = ppval(A_HR_pchip,HR/60);
B_HR_pchip = pchip(freq_all,para_fitted_Ca(3,:));
B_HR = ppval(B_HR_pchip,HR/60);
C_HR_pchip = pchip(freq_all,para_fitted_Ca(4,:));
C_HR = ppval(C_HR_pchip,HR/60);
Ca0_HR_pchip = pchip(freq_all,para_fitted_Ca(5,:));
Ca0_HR = ppval(Ca0_HR_pchip,HR/60);

stim_period = 1/(HR/60);

xm_LV   = -0.60;
xm_SEP  = 0.40;
xm_RV   = 1.0;
ym    = 0.50;

SL_LV  = 2.2;
SL_SEP = 2.2;
SL_RV  = 2.2;

V_SA = 3.0;
V_SV = 4.80;
V_PA = 0.5;
V_PV = 1.0; 
V_Ao = 1.0;

P1_0_LV = 0; % 0th moment state A1, LV
P1_1_LV = 0; % 1st moment state A1, LV
P1_2_LV = 0; % 2nd moment state A1, LV
P2_0_LV = 0; % 0th moment state A2, LV
P2_1_LV = 0; % 1st moment state A2, LV
P2_2_LV = 0; % 2nd moment state A2, LV
P3_0_LV = 0; % 0th moment state A3, LV
P3_1_LV = 0; % 1st moment state A3, LV
P3_2_LV = 0; % 2nd moment state A3, LV
N_LV = 1;
U_NR_LV = 0;
P1_0_RV = 0; % 0th moment state A1, LV
P1_1_RV = 0; % 1st moment state A1, LV
P1_2_RV = 0; % 2nd moment state A1, LV
P2_0_RV = 0; % 0th moment state A2, LV
P2_1_RV = 0; % 1st moment state A2, LV
P2_2_RV = 0; % 2nd moment state A2, LV
P3_0_RV = 0; % 0th moment state A3, LV
P3_1_RV = 0; % 1st moment state A3, LV
P3_2_RV = 0; % 2nd moment state A3, LV
N_RV = 1;
U_NR_RV = 0;
P1_0_SEP = 0; % 0th moment state A1, LV
P1_1_SEP = 0; % 1st moment state A1, LV
P1_2_SEP = 0; % 2nd moment state A1, LV
P2_0_SEP = 0; % 0th moment state A2, LV
P2_1_SEP= 0; % 1st moment state A2, LV
P2_2_SEP = 0; % 2nd moment state A2, LV
P3_0_SEP = 0; % 0th moment state A3, LV
P3_1_SEP = 0; % 1st moment state A3, LV
P3_2_SEP = 0; % 2nd moment state A3, LV
N_SEP = 1;
U_NR_SEP = 0;

init = [xm_LV ,xm_SEP ,xm_RV ,ym , SL_LV, SL_SEP, SL_RV, V_LV, V_RV, ...
       P1_0_LV, P1_1_LV, P1_2_LV ,P2_0_LV, P2_1_LV, P2_2_LV, P3_0_LV, P3_1_LV, P3_2_LV, N_LV, U_NR_LV,...
       P1_0_SEP,P1_1_SEP,P1_2_SEP,P2_0_SEP,P2_1_SEP,P2_2_SEP,P3_0_SEP,P3_1_SEP,P3_2_SEP,N_SEP,U_NR_SEP,...
       P1_0_RV, P1_1_RV, P1_2_RV, P2_0_RV, P2_1_RV, P2_2_RV, P3_0_RV, P3_1_RV, P3_2_RV, N_RV, U_NR_RV,...
       V_SV, V_PV ,V_SA ,V_PA, V_Ao]';

opts = optimset('Display','iter','MaxFunEvals',100000,'MaxIter',10000);
TrisegEquations(init(1:4),Vw_LV,Vw_SEP,Vw_RV,SL_LV,SL_SEP,SL_RV,V_LV,V_RV,Amref_LV,Amref_SEP,Amref_RV);
x = fsolve(@TrisegEquations,init(1:4),opts,Vw_LV,Vw_SEP,Vw_RV,SL_LV,SL_SEP,SL_RV,V_LV,V_RV,Amref_LV,Amref_SEP,Amref_RV);
init(1:4) = x;
% Lumped circulatory parameters
C_Ao = 0.0022045;  % Proximal aortic compliance, mL/mmHg
C_SA = 0.0077157; % Systemic arterial compliance, mL/mmHg
C_SV = 2.5; % Systemic venous compliance, mL/mmHg  DAB 10/7/2018
C_PV = 0.25; % Pulmonary venous compliance, mL/mmHg
C_PA = 0.013778; % Pulmonary arterial compliance, mL/mmHg
R_Ao   = 2.5; % resistance of aorta , mmHg*sec/mL
R_SA   = adjvar(7)*88/CO_target*60;% mmHg*sec/mL; % Systemic vasculature resistance, mmHg*sec/mL
% R_SA   = 2.25*88/CO_target*60;% mmHg*sec/mL; %  TAC #1
R_PA   = 12/CO_target*60; % Pulmonary vasculature resistance, mmHg*sec/mL % Match the old code(9/5 BM) DAB change 9/15
R_SV   = 0.25; 
R_PV   = 0.25; 
R_vlv  = 0.05; %  valve resistance, mmHg*sec/mL
R_AV   = R_vlv + R_TAC; % resistance across aortic valve
R_tAo  = 0.5;
R_tSA  = 4;
Kse    = 50000; % series element elastance, mmHg/micron (Changed to match the value in Tewari's code) (9/5 BM)

M = speye(47);
M(1,1) = 0;
M(2,2) = 0;
M(3,3) = 0;
M(4,4) = 0; 
input = [CO_target stim_period Vw_LV Vw_SEP Vw_RV R_TAC MgATP_LV MgADP_LV Pi_LV MgATP_SEP MgADP_SEP Pi_SEP MgATP_RV MgADP_RV Pi_RV A_HR B_HR C_HR Ca0_HR Amref_LV Amref_SEP Amref_RV];
options = odeset('Mass',M,'RelTol',1e-6,'AbsTol',1e-6,'MaxStep',stim_period/50);
[ts,ys] = ode15s(@dXdT_cardiovascular_mechanics,[0 120*stim_period],init,options,adjvar,input);
% V_LV   = ys(:,8); % volume LV, mL
% V_RV   = ys(:,9); % volume RV, mL
% V_SV   = ys(:,43); % volume of systemic veins
% V_PV   = ys(:,44); % volume of pulmonary veins
% V_SA   = ys(:,45); % volume of systemic arterys
% V_PA   = ys(:,46); % volume of pulmonary arterys
% V_Ao   = ys(:,47); % volume of proximal aorta
% V_T = V_LV + V_RV + V_SV + V_PV + V_SA + V_PA + V_Ao;
% 
% P_PV = V_PV/C_PV;
% P_PA = V_PA/C_PA;
% 
%  figure(1); plot(ts,V_LV,ts,ys(:,11)); title('ventricular volumes')
% % figure(2); plot(ts,(V_PV + V_PA)./V_T);
% figure(3); plot(ts,P_PA,ts,P_PV); title('pulmonary pressures');
init = ys(end,:);
[t,Y] = ode15s(@dXdT_cardiovascular_mechanics,[0 1*stim_period],init,options,adjvar,input);

% % Assignig the solution of the ODE's to the variables

xm_LV  = Y(:,1); % LV heart geometry variable, cm
xm_SEP = Y(:,2); % septum heart geometry variable, cm
xm_RV  = Y(:,3); % RV heart geometry variable, cm
ym     = Y(:,4); % Heart geometry variable, cm
SL_LV  = Y(:,5); % sarcomere length, LV, micron
SL_SEP = Y(:,6); % sarcomere length, septum, micron
SL_RV  = Y(:,7); % sarcomere length, RV, micron
V_LV   = Y(:,8); % volume LV, mL
V_RV   = Y(:,9); % volume RV, mL
% 
% P1_0_LV = Y(:,10); % 0th moment state A1, LV
% P1_1_LV = Y(:,11); % 1st moment state A1, LV
% P1_2_LV = Y(:,12); % 2nd moment state A1, LV
% P2_0_LV = Y(:,13); % 0th moment state A2, LV
% P2_1_LV = Y(:,14); % 1st moment state A2, LV
% P2_2_LV = Y(:,15); % 2nd moment state A2, LV
P3_0_LV = Y(:,16); % 0th moment state A3, LV
P3_1_LV = Y(:,17); % 1st moment state A3, LV
P3_2_LV = Y(:,18); % 2nd moment state A3, LV
N_LV    = Y(:,19); % non-permissive fraction LV
U_NR_LV = Y(:,20); % U_NR represents the Non relaxed state

% P1_0_SEP = Y(:,21); % 0th moment state A1, SEP
% P1_1_SEP = Y(:,22); % 1st moment state A1, SEP
% P1_2_SEP = Y(:,23); % 2nd moment state A1, SEP
% P2_0_SEP = Y(:,24); % 0th moment state A2, SEP
% P2_1_SEP = Y(:,25); % 1st moment state A2, SEP
% P2_2_SEP = Y(:,26); % 2nd moment state A2, SEP
P3_0_SEP = Y(:,27); % 0th moment state A3, SEP
P3_1_SEP = Y(:,28); % 1st moment state A3, SEP
P3_2_SEP = Y(:,29); % 2nd moment state A3, SEP
N_SEP    = Y(:,30); % nonpermissive fraction SEP
U_NR_SEP = Y(:,31); % U_NR represents the Non relaxed state

% P1_0_RV = Y(:,32); % 0th moment state A1, RV
% P1_1_RV = Y(:,33); % 1st moment state A1, RV
% P1_2_RV = Y(:,34); % 2nd moment state A1, RV
% P2_0_RV = Y(:,35); % 0th moment state A2, RV
% P2_1_RV = Y(:,36); % 1st moment state A2, RV
% P2_2_RV = Y(:,37); % 2nd moment state A2, RV
P3_0_RV = Y(:,38); % 0th moment state A3, RV
P3_1_RV = Y(:,39); % 1st moment state A3, RV
P3_2_RV = Y(:,40); % 2nd moment state A3, RV
N_RV    = Y(:,41); % nonpermissive fraction RV
U_NR_RV = Y(:,42); % U_NR represents the Non relaxed state

V_SV   = Y(:,43); % volume of systemic veins
V_PV   = Y(:,44); % volume of pulmonary veins
V_SA   = Y(:,45); % volume of systemic arterys
V_PA   = Y(:,46); % volume of pulmonary arterys
V_Ao   = Y(:,47); % volume of proximal aorta
V_T = V_LV + V_RV + V_SV + V_PV + V_SA + V_PA + V_Ao;

% PlotTriSeg(xm_LV,xm_SEP,xm_RV,ym,t)
%  Pulmonary Pressures
P_PV = V_PV/C_PV;
P_SV = V_SV/C_SV;
P_PA = V_PA/C_PA;
P_SA = V_SA/C_SA;

Am_LV = pi*(xm_LV.^2 + ym.^2);
Am_SEP = pi*(xm_SEP.^2 + ym.^2);
Am_RV = pi*(xm_RV.^2 + ym.^2);
Cm_LV = 2*xm_LV./(xm_LV.^2 + ym.^2);
Cm_SEP = 2*xm_SEP./(xm_SEP.^2 + ym.^2);
Cm_RV = 2*xm_RV./(xm_RV.^2 + ym.^2);
z_LV = 3*Cm_LV.*Vw_LV./(2*Am_LV);
z_SEP = 3*Cm_SEP.*Vw_SEP./(2*Am_SEP);
z_RV = 3*Cm_RV.*Vw_RV./(2*Am_RV);

epsf_LV = (1/2)*log(Am_LV./Amref_LV) - (1/12)*z_LV.^2 - 0.019*z_LV.^4;
epsf_SEP = (1/2)*log(Am_SEP./Amref_SEP) - (1/12)*z_SEP.^2 - 0.019*z_SEP.^4;
epsf_RV = (1/2)*log(Am_RV./Amref_RV) - (1/12)*z_RV.^2 - 0.019*z_RV.^4;
SLo_LV = Lsref*exp(epsf_LV);
SLo_SEP = Lsref*exp(epsf_SEP);
SLo_RV = Lsref*exp(epsf_RV);

% % Total forces
sigmaf_LV = -Kse*(SL_LV - SLo_LV);
sigmaf_SEP = -Kse*(SL_SEP - SLo_SEP);
sigmaf_RV = -Kse*(SL_RV - SLo_RV);

% % equilibrium of forces at junction circle
Tm_LV = (Vw_LV.*sigmaf_LV./(2*Am_LV)).*(1 + z_LV.^2/3 + z_LV.^4/5);
% Tm_SEP = (Vw_SEP.*sigmaf_SEP./(2*Am_SEP)).*(1 + z_SEP.^2/3 + z_SEP.^4/5);
Tm_RV = (Vw_RV.*sigmaf_RV./(2*Am_RV)).*(1 + z_RV.^2/3 + z_RV.^4/5);
sinalpha_LV = 2*xm_LV.*ym./(xm_LV.^2 + ym.^2);
% sinalpha_SEP = 2*xm_SEP.*ym./(xm_SEP.^2 + ym.^2);
sinalpha_RV = 2*xm_RV.*ym./(xm_RV.^2 + ym.^2);
% cosalpha_LV = (-xm_LV.^2 + ym.^2)./(xm_LV.^2 + ym.^2);
% cosalpha_SEP = (-xm_SEP.^2 + ym.^2)./(xm_SEP.^2 + ym.^2);
% cosalpha_RV = (-xm_RV.^2 + ym.^2)./(xm_RV.^2 + ym.^2);
Tx_LV = Tm_LV.*sinalpha_LV;
% Tx_SEP = Tm_SEP.*sinalpha_SEP;
Tx_RV = Tm_RV.*sinalpha_RV;
% Ty_LV = Tm_LV.*cosalpha_LV;
% Ty_SEP = Tm_SEP.*cosalpha_SEP;
% Ty_RV= Tm_RV.*cosalpha_RV;

% % ventricular pressure
ptrans_LV = 2*Tx_LV./ym;
ptrans_RV = 2*Tx_RV./ym;
P_LV = -ptrans_LV;
P_RV = ptrans_RV;
% Ao valves closed equations
P_Ao_closed = (C_SA*R_Ao*R_SA*V_Ao + C_SA*R_Ao*R_tSA*V_Ao + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));
P_SA_closed = (C_Ao*R_Ao*R_SA*V_SA + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA + C_Ao*C_SA*P_SV*R_Ao*R_tSA + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));
% Ao valve open equations 
P_Ao_open = (C_SA*R_Ao*R_SA*R_AV*V_Ao + C_SA*R_Ao*R_tSA*R_AV*V_Ao + C_SA*R_SA*R_tSA*R_AV*V_Ao + C_Ao*R_SA*R_tAo*R_AV*V_SA + C_Ao*C_SA*P_LV*R_Ao*R_SA*R_tAo + C_Ao*C_SA*P_LV*R_Ao*R_tSA*R_tAo + C_Ao*C_SA*P_LV*R_SA*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo*R_AV)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_AV + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_AV + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_AV + R_SA*R_tAo*R_AV + R_tSA*R_tAo*R_AV));
P_SA_open = (C_Ao*R_Ao*R_SA*R_tAo*V_SA + C_Ao*R_Ao*R_SA*R_AV*V_SA + C_SA*R_SA*R_tSA*R_AV*V_Ao + C_Ao*R_SA*R_tAo*R_AV*V_SA + C_Ao*C_SA*P_SV*R_Ao*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_Ao*R_tSA*R_AV + C_Ao*C_SA*P_LV*R_SA*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo*R_AV)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_AV + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_AV + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_AV + R_SA*R_tAo*R_AV + R_tSA*R_tAo*R_AV));
QOUT_LV = -(C_SA*R_Ao*R_SA*V_Ao + C_SA*R_Ao*R_tSA*V_Ao + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA - C_Ao*C_SA*P_LV*R_Ao*R_SA - C_Ao*C_SA*P_LV*R_Ao*R_tSA - C_Ao*C_SA*P_LV*R_SA*R_tSA - C_Ao*C_SA*P_LV*R_SA*R_tAo - C_Ao*C_SA*P_LV*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_AV + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_AV + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_AV + R_SA*R_tAo*R_AV + R_tSA*R_tAo*R_AV));

P_Ao = P_Ao_open.*(P_LV>P_Ao_open) + ...
       P_Ao_closed.*(P_LV<=P_Ao_open);
QOUT_LV = QOUT_LV.*(P_LV>P_Ao_open);

CO_sim =  (max(V_LV)-min(V_LV))*HR;

CO_mean = mean(QOUT_LV);

SV_LV_sim = max(1e3*V_LV) - min(1e3*V_LV);
EF_LV_sim = SV_LV_sim/max(1e3*V_LV) * 100;
SV_RV_sim = max(1e3*V_RV) - min(1e3*V_RV);
EF_RV_sim = SV_RV_sim/max(1e3*V_RV) * 100;
    
%% calculate the error
edLV_sim =  max(1e3*V_LV);
esLV_sim =  min(1e3*V_LV);
edRV_sim =  max(1e3*V_RV);
esRV_sim =  min(1e3*V_RV);

g2_LV = (MgATP_LV/K_T)/(1.0 + MgATP_LV/K_T + MgADP_LV/K_D);
k3_LV = k3*g2_LV;%
g2_SEP = (MgATP_SEP/K_T)/(1.0 + MgATP_SEP/K_T + MgADP_SEP/K_D);
k3_SEP = k3*g2_SEP;%
g2_RV = (MgATP_RV/K_T)/(1.0 + MgATP_RV/K_T + MgADP_RV/K_D);
k3_RV = k3*g2_RV;%
f_alpha3o_LV  = (P3_0_LV + alpha3*(s3*s3*P3_0_LV + 2.0*s3*P3_1_LV + P3_2_LV)); 
f_alpha3o_SEP = (P3_0_SEP + alpha3*(s3*s3*P3_0_SEP + 2.0*s3*P3_1_SEP + P3_2_SEP)); 
f_alpha3o_RV  = (P3_0_RV + alpha3*(s3*s3*P3_0_RV + 2.0*s3*P3_1_RV + P3_2_RV)); 

% detachment rates
ti = 0:0.00001:stim_period;
MAP = mean(interp1(t,P_Ao,ti));

r_LV  = interp1(t,k3_LV*f_alpha3o_LV,ti);
r_SEP = interp1(t,k3_SEP*f_alpha3o_SEP,ti);
r_RV  = interp1(t,k3_RV*f_alpha3o_RV,ti);

% LV X-bridge turnover rate
Vw_LV_W = (2/3)*LVW/1000;
Vw_SEP_W= (1/3)*LVW/1000;
rate_of_XB_turnover_ave = (Vw_LV_W*mean(r_LV) + Vw_SEP_W*mean(r_SEP))/(Vw_LV_W + Vw_SEP_W) 

% unit convert to oxygen consumption
ATP_ase_mechannics_Averge_LV_SEP = (1.327/5.1253)*rate_of_XB_turnover_ave %  1.31 Kstiff - ATP hydrolized (mmol/s/(L cell)) per X-bridge turnover rate in LV

Fitting_error(1) = (edLV_target - max(1e3*V_LV))^2 / (edLV_target * max(1e3*V_LV)); % check error on edLV
Fitting_error(2) = ((esLV_target - min(1e3*V_LV))^2 / (esLV_target * min(1e3*V_LV)));% check error on esLV
Fitting_error(3) = ( MAP_target- MAP)^2 / (MAP_target * MAP); % MAP error
Fitting_error(4) = ((SV_LV_target - SV_LV_sim)^2 / (SV_LV_target * SV_LV_sim)); % check relative error on stroke volume

% [MAP max(V_LV) min(V_LV) P_LV(end)]

if flag_plot_figure == 1
%% Plotting
h1 = figure(1);
subplot(2,3,1); plot(t,V_LV,t,V_RV); title('ventricular volumes')
ylabel('$V$ (ml)','interpreter','latex','fontsize',18)
xlabel('$t$ (sec)','interpreter','latex','fontsize',18)

% figure(2); plot(t,(V_PV+V_PA)./V_T); title('pulmonary volume fraction');
subplot(2,3,2); plot(t,P_RV,t,P_PA,t,P_SV); title('pulmonary pressures'); legend('P_{RV}','P_{PA}','P_{SV}');  set(gca,'Ylim',[0 55]);
ylabel('$P$ (mmHg)','interpreter','latex','fontsize',18)
xlabel('$t$ (sec)','interpreter','latex','fontsize',18)

subplot(2,3,3); plot(t,Y(:,5:7)); title('SL'); legend('LV','SEP','RV'); % SL's
ylabel('$SL$ (um)','interpreter','latex','fontsize',18)
xlabel('$t$ (sec)','interpreter','latex','fontsize',18)


subplot(2,3,4); plot(t,P_LV,t,P_Ao,t,P_SA,t,P_PV);
hold on; plot([0 1/(HR/60)],[max(P_Ao) max(P_Ao)],'k--',[0 1/(HR/60)],[max(P_Ao)+4*delta_p max(P_Ao)+4*delta_p],'k--'); hold off
legend('P_{LV}','P_{Ao}','P_{SA}','P_{PV}','Max P_{Ao}','Max P_{LV}'); %set(gca,'Ylim',[0 125]);
ylabel('$P$ (mmHg)','interpreter','latex','fontsize',18)
xlabel('$t$ (sec)','interpreter','latex','fontsize',18)

subplot(2,3,5); plot(V_LV,P_LV,V_RV,P_RV,[edLV_target edLV_target]./1000,[0 130],'k--',[esLV_target esLV_target]./1000,[0 130],'k--'); set(gca,'Xlim',[0 0.75]);
ylabel('$P$ (mmHg)','interpreter','latex','fontsize',18)
xlabel('$V$ (ml)','interpreter','latex','fontsize',18)

txt_EF = ['EF =', num2str(EF_LV_sim)]
subplot(2,3,6);cla; text(0.2,0.8,txt_EF)
txt_SV = ['SV =', num2str(SV_LV_sim),' ml'];
text(0.2,0.9,txt_SV)
title(['MAP = ', num2str(MAP),' mmHg']);
set(h1,'Position',[50 50 1250 600])


if flag_swap_metabolite == 1
         filename = strcat(pwd,'\results\rest\rat_swap',num2str(rat_number));
    else 
         filename = strcat(pwd,'\results\rest\rat',num2str(rat_number)); 
end
save(filename)
end
% save data to excel file 
result_table(rat_number,:) =[rat_number, Pi_LV, MgATP_LV, MgADP_LV, EF_LV_sim, MAP, P_LV(end), dGrATPase_LV, CO_sim]
toc
end

if flag_swap_metabolite == 1
        excel_file_name = 'Result_table_swap.xlsx';
        xlRange = 'A2:I21';
        xlswrite(excel_file_name,result_table,xlRange)
    else 
        excel_file_name = 'Result_table_rest.xlsx';
        xlRange = 'A2:I21';
        xlswrite(excel_file_name,result_table,xlRange)
end