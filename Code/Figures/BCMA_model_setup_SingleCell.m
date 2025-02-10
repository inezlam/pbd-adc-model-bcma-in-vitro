%{ 
Description: Parameter and initial conditions setup for anti-BCMA-PBD ADC model.
Author: Inez Lam

Model Units = change in concentration per unit time (nM/hr)
- Model Output Units: Concentration (nM = nmol/L)
- Time Units: Hours (hr)

Inputs:
- TIME = time (optional, may be omitted by indicating ~)
- DOSETYPE = type of dose administered

Note:
- Any parameter going into structure p must be able to be "varied" in a
sensitivity analysis.

This version includes the following added mechanisms over the base model:
- Naked antibody
- Soluble antigen
- Warhead tracking
- Cleared ADC and Ab
- DNA crosslinking
%}

run('DefineParameterStuctOrder.m')

NA = 6.02e23; % (molecules/mole) Avogadro's number

%% Experimental Volumes and Number of Cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hartley 2018 - Experimental conditions for DNA Crosslinking experiments
V = 100e-6; % Volume per well in L for cytotoxicity assay Hartley 2018 paper
cellDensity = 1e8; % 1e5 cells/mL * (1e3 mL/L) = 1e8 cells per L 
Ncell = cellDensity * V;

%% Cell Volumes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vcell = 1e-12; % Volume of a cell
NucleusFraction = 0.1; % nucleus fraction of total cell volume
EndoLysoFraction = 0.05; % endosome/lysosome fraction of total cell volume

% Total Subcellular Compartment Volumes
p.Vmedia = V - (Vcell * Ncell); % total volume of media not including cells
p.Vendo = (Vcell * EndoLysoFraction); % volume of endosomes
p.Vnuc = (Vcell * NucleusFraction); % total volume of nucleus across all cells
p.Vcytoplasm = Vcell * (1 - NucleusFraction - EndoLysoFraction);

%% %%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Infusion of ADC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.q = 0; % (nmol/hr) infusion rate of ADC, assume no infusion for now

% Drug to Antibody Ratio %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.DAR = 2;

% ADC Binding to Antigen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bound Ag => Apparent Kd = 3.1 nM;
p.konAg = 1.34; % (nM*hr)^-1 binding rate constant of ADC to Ag
p.koffAg = 4.16268;  % 1/hr unbinding rate constant of ADC-Ag Complex

% Soluble Ag => Kd = 60.7;
p.konAg_s = 1.34; % (nM*hr)^-1 - measured in ONC2228-0002 (monomeric BCMA)
p.koffAg_s = 81.72;  % 1/hr - measured in ONC2228-0002 (monomeric BCMA)

% Internalization and Recycling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Of Complex
p.ki = 0.094367;
p.frec = 0.43444;
p.ke = 2.9612;

%% Of Antigen
p.kiAg = p.ki;
p.frecAg = p.frec;
p.keAg = p.ke;
p.kAgProduction = 1.67842e-12;

%% Of ADC
p.kiADC = 1.2095;
p.keADC = 25;
p.krecADC = 9.2597;

%% FcRn Binding
p.koffFcRn = p.koffAg_s; % 1/hr
p.konFcRn = p.koffFcRn/648; % (nM*hr)^-1 % Calculated from known KD of MEDI2228 to FcRn

%% Antigen Production <--CHECK
% Degradation of ADC Extracellularly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.kdeconjugation = 0; % (hr-1) ADC lysosomal degradation rate constant
kdeconjugation_MouseSerum = 3.2617e-04; % (hr-1) ADC lysosomal degradation rate constant

% Efflux and Influx %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.kinf = log(2)/3.6 * 60; % p.kinf = 11.55; % (hr-1) influx rate constant
R = 82.45;
p.keff = p.kinf/(1 + R); % p.keff = 0.1384; % (hr-1) efflux rate constant

% Trafficking to Nucleus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.kinN = 0.1384; %0.092; % rate constant for trafficking into the nucleus
p.koutN = 11.55; %0.092; % rate constant for trafficking out of the nucleus

% Warhead Binding to DNA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.konDNA = 0.06; %1e4 * 60; % (nM^-1 hr^-1) binding rate constant of Wn to DNA
p.koffDNA = 0.006; %1e-5 * 60; % (1/hr) unbinding rate constant of Wn-DNA complex

% DNA Crosslinking (Hill Coefficient) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.n = 1; % Hill coefficient
Crosslink50_Dose = 1.2398; % (nM) PBD dose causing 50% crosslinking 
p.KA = 121.6738; % (nM) Wn-DNA concentration causing 50% crosslinking (calculated)

% Cell Growth and Inhibition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
doublingTime = 50; % doubling time for NCI-H929
p.kg = log(2)/doublingTime; % (hr-1) 

% Fitted Parameter
p.kkill = 0.11708; % (hr-1) warhead cell killing rate constant for SG3199

%% %%%%%%%%%%%%%%%%%%%%%%%%% Initial Amounts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ADC_MW = 151.614 * (10^3); % kDa = kg/mol => g/mol % MW of MEDI2228
ADC_doses = [1.85999E-05 5.58656E-05 0.000167531 0.000502644 0.001507935 0.004523797 0.013571392 0.040714182 0.122142546 0.366427639 1.099283048 3.297848484]; % nM - initial concentrations of ADC
ADCDoseConversionFactor = (1e-3) * (1e9) / ADC_MW; % ug/mL => nmol/L (nM)
ADC_doses_9 = [7.62079E-05 0.000228624 0.000685871 0.002057613 0.00617284 0.01851852 0.05555556 0.1666667 0.5].*ADCDoseConversionFactor; % ug/mL => nmol/L (nM)
Ag0 = (18931) / NA * 10^9 ; % % Ag on cell surface - volume of fluid in well (molecules/cell => nmol/cell)
Ag_endo0 = 20.04285 * p.Vendo; % calculated Ag in endosome (nmol)
DNA0 = (3e9/10)/2*1 / NA * 10^9 * 0.01; % initial DNA (nmol)
Ag_s_MW = 5300; % soluble Ag molecular weight (g/mol)
Ag_s0 = 0; % initial soluble Ag (nmol)
FcRn0 = 4.98e4 * p.Vendo; % initial FcRn (nmol)

%% %%%%%%%%%%%%%%%%%%%%%% Initial Conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numVars = 36;
eqns = num2cell(1:numVars);
[ADC, Ag, ADC_Ag, Ag_s, ADC_Ag_s,...
    Ab, Ab_Ag, Ab_Ag_s,...
    ADC_endo, Ag_endo, ADC_Ag_endo,...
    FcRn_endo,ADC_FcRn_endo,...
    Ab_endo, Ab_Ag_endo, Ab_FcRn_endo,...
    ADC_lys, Ag_lys, ADC_Ag_lys,...
    Ab_lys, Ab_Ag_lys,...
    ADC_cl, Ag_cl, Ab_cl,...
    Wex, Wex_deg, Wex_efflux, Win, Win_lys, Win_deg, Win_influx, Win_nuc,...
    Wn, DNA, Wn_DNA, Cells] = deal(eqns{:});
varNames = {'ADC','Ag','ADC-Ag','Ag_s','ADC-Ag_s',...
    'Ab','Ab-Ag','Ab-Ag_s',...
    'ADC_e_n_d_o','Ag_e_n_d_o','ADC-Ag_e_n_d_o',...
    'FcRn_e_n_d_o','ADC-FcRn_e_n_d_o',...
    'Ab_e_n_d_o','Ab-Ag_e_n_d_o','Ab-FcRn_e_n_d_o',...
    'ADC_l_y_s','Ag_l_y_s','ADC-Ag_l_y_s',...
    'Ab_l_y_s','Ab-Ag_l_y_s',...
    'ADC_c_l','Ag_c_l','Ab_c_l',...
    'W_e_x','W_e_xdeg','W_e_xefflux','W_i_n','W_i_nlys','W_i_ndeg','W_i_ninflux','W_i_nnuc',...
    'Wn','DNA','Wn-DNA','Cells'};

%% Set initial conditions
InitCond = zeros(numVars,1);
InitCond(Ag) = Ag0/p.Vmedia; % Initial extracellular membrane-bound Ag 
InitCond(Ag_s) = Ag_s0/p.Vmedia; % Initial extracellular soluble Ag
InitCond(Ag_endo) = Ag_endo0/p.Vendo; % Initial intracellular Ag 
InitCond(Ag_lys) = 33.566487;
InitCond(FcRn_endo) = FcRn0/p.Vendo;
InitCond(DNA) = DNA0/p.Vnuc;
InitCond(Cells) = Ncell; % Concentration of cells in nM

InitCond_noADC = InitCond; % nM
InitCond_noADC(ADC) = 0;

InitCond_PBD = InitCond_noADC;
InitCond_PBD(Wex) = 1.69;

% Set conditions depending on dose type
if doseType == "ADC"
    InitCond(ADC) = ADC_doses(8);
elseif doseType == "Ab"
    InitCond(Ab) = ADC_doses(8);
elseif doseType == "PBD"
    InitCond(ADC) = 0;
    InitCond(Wex) = 1.69;
elseif doseType == "IsotypeADC"
    InitCond(ADC) = ADC_doses(8);
    % Isotype ADC
    p.konAg = 0; % (nM*hr)^-1 binding rate constant of ADC to Ag
    p.koffAg = 0;  % 1/hr unbinding rate constant of ADC-Ag Complex
    p.konAg_s = 0; % (nM*hr)^-1
    p.koffAg_s = 0;  % 1/hr
else 
    InitCond(ADC) = 0; 
    InitCond(Ab) = 0;
    InitCond(Wex) = 0;
end

%% Optimize for kAgProduction and Calculate Ag_endo_ss
calculateAgProduction = "No"; % "Yes" or "No"
plotAgProdFigure = "No"; % "Yes" or "No" 
calculateAg_endo_ss_only = "Yes"; % "Yes" or "No" 
if calculateAgProduction == "Yes"
    run("AgProd_Optimization.m");
    p.kAgProduction = optimal_AgProduction; % Set new kAgProduction value
    % Set Ag_endo_ss as initial condition for Ag_endo
    InitCond(Ag_endo) = Y0(end,Ag_endo);
    InitCond_noADC(Ag_endo) = Y0(end,Ag_endo);
    InitCond_PBD(Ag_endo) = Y0(end,Ag_endo);
    InitCond_noADC(Ag_lys) = Y0(end,Ag_lys);
    InitCond(Ag_lys) = Y0(end,Ag_lys);
    InitCond_PBD(Ag_lys) = Y0(end,Ag_lys);

elseif calculateAg_endo_ss_only == "Yes"
    fprintf("Running pre-simulation to determine steady state Ag values...")   
    % Calculate Ag_endo_ss only
    time_ss = 24*28;
    [T_ss,Ag_ss_calculated,Ag_endo_ss,Ag_lys_ss] = Ag_ss_calculation(eqns_file,ode_options,p,time_ss,Ag,Ag_endo,Ag_lys,InitCond_noADC,p.kAgProduction);
    % Set Ag_endo_ss as initial condition for Ag_endo
    InitCond_noADC(Ag_endo) = Ag_endo_ss(end);
    InitCond(Ag_endo) = Ag_endo_ss(end);
    InitCond_PBD(Ag_endo) = Ag_endo_ss(end);
    InitCond_noADC(Ag_lys) = Ag_lys_ss(end);
    InitCond(Ag_lys) = Ag_lys_ss(end);
    InitCond_PBD(Ag_lys) = Ag_lys_ss(end);

    % % Plot comparison
    % figure;
    % set(gcf,'color','w')
    % yline(InitCond(Ag),'r--','linewidth',2)
    % hold on;
    % plot(T_ss,Ag_ss_calculated,'linewidth',2)
    % title({'k_{AgProduction} at Steady State = ',[num2str(p.kAgProduction) ' nmol/hr']},'fontsize',20)
    % ylabel('Concentration (nM)','fontsize',20)
    % xlabel('Time (hours)','fontsize',20)
    % legend({'Ag_s_s','Ag'})
    % set(gca,'FontSize',20)
    % grid on; box on;
    fprintf("Ag_endo_ss is %e nM, Ag_lys_ss is %e nM.\n",Ag_endo_ss(end),Ag_lys_ss(end)) 
end  

%% Calculate KA value
calculateKA = "Yes"; % "Yes" or "No"
plotCrosslinkFigure = "No"; % "Yes" or "No" 
if calculateKA == "Yes"
    fprintf("Calculating value of KA (half maximal concentration of Wn-DNA for crosslinking)... ")    
    ExposureTime = 2;
    IncubationTime = 24;
    [p.KA,T_KA,Y_KA] = KA_calculation(eqns_file,p,ExposureTime,IncubationTime,Crosslink50_Dose,InitCond_noADC,Wex,Wn_DNA,plotCrosslinkFigure);
    fprintf("The calculated value for KA is %e nM.\n",p.KA)
end  
