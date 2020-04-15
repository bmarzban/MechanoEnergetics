% This code is the main driver for the cardiovascular mechanics model 
clear; 
tic;
%% chooosing the rat number (Mean Sham rat is Rat number 9,Mean TAC rat is number 19 )
rat_number = 19

if rat_number<=9
    shamRat = 1;
else
    shamRat = 0;
end

%% Read the experimental data for SHAM and TAC rats from the excel file 
data = xlsread('data1.xlsx','A3:W23');
BW  = data(rat_number , 1); % g
LVW = data(rat_number , 2); % mg
RVW = data(rat_number , 3); % mg
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
if rat_number<=9
%     shamRat = 1;
Ox_capacity = Ox_capacity_sham;
end

% Average sham
TAN_sham = data(9,16)/1000; % mole/L cell
CRtot_sham = data(9,18)/1000; % mole/L cell
% TEP_sham =1.6* data(9,20)/1000; % mole/L cell

A0 = 8.62; %mmol/ (L Cell)
P0 = 29.78;%mmol/ (L Cell)
a = 0.082; %mmol/ (L Cell. year)
p = 0.283; %mmol/ (L Cell. year)

TEP = (P0 - p*((A0 - TAN*1000) /a))/1000; % mole/L cell
TEP_sham = (P0 - p*((A0 - TAN_sham*1000) /a))/1000; %mole/L cell

if shamRat == 0
    preV = data(rat_number , 11); % mm/s
    postV = data(rat_number , 12); % mm/s
    postV = postV/1000; %m/s
    preV = preV/1000; %m/s
    rho_blood = 1060; % kg/m^3
    delta_p = 0.5*(postV^2-preV^2)*rho_blood; % Pa
    delta_p = 0.0075*delta_p; % mmHg
    if rat_number == 19
        delta_p = 31.48;
    end
    R_TAC = delta_p/CO_target*60;
else 
    R_TAC = 0;
end

CO_target = 95; %change the Co to 95 to have the same Resitance paramters. as mean sham , we had that the CO changed since we wanted to have the 
% right R_TAC
%% Adjustable variables
%(Tune the following variables to fit the EDLV, ESLV, EDRV, ESLV, CO, ATP consumption Rate predicted by Ox-Phos and Mechanic model)

%         and   x_ATPase ==> ATP hydrolsis rate for LV, Septal, and RV independantly  

% adjvar = [Reference area LV & Septal,  Reference area RV,  k_passive, ksr, kforce, R_SA, R_TAC]

%% para set 4
 
adjvar = [1.105 0.93 1.34 2.0812*1.11  1.1*0.81 1.1*0.81 2.28 0.54]; % Rat 19 1.31 kstiff % eta = 0.1(19 is mean TAC rat)
tune_ATPase_LV = 0.861* (1/ 0.6801) *1.0e-3;

R_TAC = adjvar(8)*R_TAC;

tune_ATPase_SEP = tune_ATPase_LV;
tune_ATPase_RV =  tune_ATPase_LV;

Amref_LV  = adjvar(1) * 2.077 ; % LV midwall reference surface area, cm^2
Amref_SEP = adjvar(2) * Amref_LV * 0.590 ; % SEP midwall reference surface area, cm^2
Amref_RV  = adjvar(3) * 3.3 ; % RV midwall reference surface area, cm^2

Vw_LV = (LVW*2/3)/1000/1.05;
Vw_SEP =(LVW/3)/1000/1.05;
Vw_RV = RVW/1000/1.05;

V_LV  = edLV_target/1000;% intial value for V_LV and V_RV assumed to be equal to edLV_target
V_RV  = edLV_target/1000;%
% V_LV  = 0.600;% intial value for V_LV and V_RV assumed to be equal to edLV_target
% V_RV  = 0.600;%
% % 
V_SA = 3.0;
V_SV = 4.80;
V_PA = 0.5;
V_PV = 1.0; 
V_Ao = 1.0;


%% parameters for calculating ATPase rate
Lsref = 1.9;
k3      = 144.5586; % transition A3 to P rate constant, 1/sec
K_T = 0.4897; 
K_D = 0.194;% Used the values from Tewari etal JMCC (9/5 BM)
alpha3  = 0.1*59.3; % Stretch sensing parameter for k3, 1/um
s3      = 9.9e-3;  % Strain at which k3 is minimum, um


%% Run the energetics model to get the metabolite concentrations

    energtics_output_LV  = EnergeticsModelScript(TAN, CRtot, TEP, Ox_capacity, tune_ATPase_LV);
%     energtics_output_SEP = EnergeticsModelScript(TAN, CRtot, TEP, Ox_capacity, tune_ATPase_SEP);
%     energtics_output_RV  = EnergeticsModelScript(TAN, CRtot, TEP, Ox_capacity, tune_ATPase_RV);
%     
%     
    MgATP_LV = energtics_output_LV(1);
    MgADP_LV = energtics_output_LV(2);
    Pi_LV = energtics_output_LV(10)*1000;

%     MgATP_SEP = energtics_output_SEP(1);
%     MgADP_SEP = energtics_output_SEP(2);
%     Pi_SEP = energtics_output_SEP(10)*1000;
%     
%     MgATP_RV = energtics_output_RV(1);
%     MgADP_RV = energtics_output_RV(2);
%     Pi_RV = energtics_output_RV(10)*1000;
% 
%    

MgATP_SEP = MgATP_LV;
MgADP_SEP = MgADP_LV;
Pi_SEP = Pi_LV ;
MgATP_RV = MgATP_LV;
MgADP_RV = MgADP_LV;
Pi_RV = Pi_LV ;

%% Run cardiovascular mechanics model

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

opts = optimset('MaxFunEvals',10000,'MaxIter',1000);
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
R_SA   =  adjvar(7)*88/CO_target*60;% mmHg*sec/mL; % Systemic vasculature resistance, mmHg*sec/mL
% R_SA   = 2.25*88/CO_target*60;% mmHg*sec/mL; %  TAC #1
R_PA   = 12/CO_target*60; % Pulmonary vasculature resistance, mmHg*sec/mL % Match the old code(9/5 BM) DAB change 9/15
R_SV   = 0.25; 
R_PV   = 0.25; 
R_vlv  = 0.05; %  valve resistance, mmHg*sec/mL
R_AV   = R_vlv + R_TAC; % resistance across aortic valve
R_tAo  = 0.5;
R_tSA  = 4;
Kse    = 50000; % series element elastance, mmHg/micron (Changed to match the value in Tewari's code) (9/5 BM)

% load init

M = speye(47);
M(1,1) = 0;
M(2,2) = 0;
M(3,3) = 0;
M(4,4) = 0; 
input = [CO_target stim_period Vw_LV Vw_SEP Vw_RV R_TAC MgATP_LV MgADP_LV Pi_LV MgATP_SEP MgADP_SEP Pi_SEP MgATP_RV MgADP_RV Pi_RV A_HR B_HR C_HR Ca0_HR Amref_LV Amref_SEP Amref_RV];
options = odeset('Mass',M,'RelTol',1e-6,'AbsTol',1e-6,'MaxStep',stim_period/50);
[t,Y] = ode15s(@dXdT_cardiovascular_mechanics,[0 60*stim_period],init,options,adjvar,input);
init = Y(end,:);
save init init
[t,Y] = ode15s(@dXdT_cardiovascular_mechanics,[0 1*stim_period],init,options,adjvar,input);
% 
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

rate_of_XB_turnover_ave = (Vw_LV*mean(r_LV) + Vw_SEP*mean(r_SEP))/(Vw_LV + Vw_SEP) 

% unit convert to oxygen consumption
% ATP_ase_mechannics_Averge_LV_SEP = (1.319/6.6079)*rate_of_XB_turnover_ave % ATP hydrolized (mmol/s/(L cell)) per X-bridge turnover rate in LV
ATP_ase_mechannics_Averge_LV_SEP = (1.327/5.1253)*rate_of_XB_turnover_ave %  1.31 Kstiff - ATP hydrolized (mmol/s/(L cell)) per X-bridge turnover rate in LV


Fitting_error(4) = (edLV_target - max(1e3*V_LV))^2 / (edLV_target * max(1e3*V_LV));
Fitting_error(5) = ((esLV_target - min(1e3*V_LV))^2 / (esLV_target * min(1e3*V_LV)));
Fitting_error(6) = (94 - MAP)^2 / (94 * MAP);
Fitting_error(8) = ((SV_LV_target - SV_LV_sim)^2 / (SV_LV_target * SV_LV_sim));
    
toc

[MAP max(V_LV) min(V_LV) P_LV(end)]
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

figure(2)
plot(t,U_NR_LV,t,U_NR_SEP,t,U_NR_RV)
title('U_NR (Non Relaxed)')
% filename = strcat(pwd,'\results\rest\rat',num2str(rat_number));
% save(filename)

% % SW = 0.0001333*polyarea(V_LV,P_LV) % Joules per beat; 1 mmHg*mL = 0.0001333 Joules
% % SW_per_min = SW * HR/ ((2/3)*LVW/1000)% % Joules per min/ per gram left ventricle. 1 mmHg*mL = 0.0001333 Joules
