%DriverGeneticAlgorithm_MechanoEnegergetic
clear; 
tic;
flag_plot_figure = 1;
flag_swap_metabolite = 0;

%% chooosing the rat number (Mean Sham rat is Rat number 9,Mean TAC rat is number 19 )
rat_number = 1;
rat_number 
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
              1.105 0.93 1.34 1.695 0.7 0.7 1 0.54 1.307 0.6 0.5; % rat 15 TAC 6- SWAP Metabolite with mean SHAM
              1.25 0.97 1.04 1.3528 0.825 0.825 1 0.5 1.366 0.6 0.5; % rat 16 TAC 7- SWAP Metabolite with mean SHAM
              1.357 0.95 1.56 1.37 0.8694 0.8694 1. 0.6 1.532 0.6 0.5; % rat 17 TAC 8- SWAP Metabolite with mean SHAM
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
SL_MAX_target = 2.2; % um

inputs = [rat_number , flag_swap_metabolite, shamRat, delta_p, LVW, RVW, HR, edLV_target, esLV_target, TAN, CRtot, Ox_capacity, TEP, R_TAC, CO_target ,MAP_target ,SL_MAX_target]

gaoptions = optimoptions('ga','MaxGenerations',50,'Display','iter');
gaoptions = optimoptions(gaoptions,'UseParallel',true);
gaoptions = optimoptions(gaoptions,'PopulationSize',100);
gaoptions = optimoptions(gaoptions,'FunctionTolerance',1e-3);
% gaoptions = optimoptions(gaoptions,'OutputFcn',@GA_DISP)
lb = adjvar - adjvar/2;
ub = adjvar + adjvar/2;

numberOfVariables =11;
ConstraintFunction = [];
ObjectiveMechanoEnergetics = @ObjectiveMechanoEnergetics;
adjvar1 = ga(ObjectiveMechanoEnergetics,numberOfVariables,[],[],[],[],lb,ub,ConstraintFunction,gaoptions)
% adjvar1 = [1.4295    0.8617    0.8780    1.7048    0.6988    0.7069    1.2274 0.6248    0.7598    0.6009    0.5305]
fitness = ObjectiveMechanoEnergetics(adjvar)
toc