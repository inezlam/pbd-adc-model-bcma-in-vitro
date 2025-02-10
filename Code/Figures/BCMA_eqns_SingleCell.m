function dydt = BCMA_eqns_SingleCell(~,y,p)
%{ 
Description: Function containing coupled, nonlinear ODE model for anti-BCMA-PBD ADCs.
Author: Inez Lam
Last Updated: 20240321

Model Units = change in concentration per unit time (nM/hr)

Inputs:
- t = time (optional, may be omitted by indicating ~)
- y = concentrations of the molecular species
- p = parameter vector

This version includes the following added mechanisms over the base model:
- Naked antibody
- Soluble antigen
- Warhead tracking
- Cleared ADC and Ab
- DNA crosslinking
%}

%% ODE Model Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Units = change in concentration per unit time (nM/hr)
dydt = zeros(length(y),1); % vector to hold ODEs

eqns = num2cell(1:length(dydt));

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

%% Equation Snippets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (+) = "taking away from", (-) = "receiving from"

% ADC-Ag: Binding/Unbindg, Internalization, Recycling, Trafficking
BindUnbind_ADC_Ag = - p.konAg*y(ADC)*y(Ag) + p.koffAg*y(ADC_Ag); % Target-Antigen Binding and Unbinding
Int_ADC_Ag = - p.ki*y(ADC_Ag); % Complex Internalization
Rec_ADC_Ag = - (p.ke*y(ADC_Ag_endo)) * p.frec; % Complex Recycling
TraffickToLys_ADC_Ag = - (p.ke*y(ADC_Ag_endo)) * (1 - p.frec); % Complex Trafficking from Endosome to Lysosome

% Ab-Ag: Binding/Unbindg, Internalization, Recycling, Trafficking
BindUnbind_Ab_Ag = - p.konAg*y(Ab)*y(Ag) + p.koffAg*y(Ab_Ag); % Target-Antigen Binding and Unbinding
Int_Ab_Ag = - p.ki*y(Ab_Ag); % Complex Internalization
Rec_Ab_Ag = - (p.ke*y(Ab_Ag_endo)) * p.frec; % Complex Recycling
TraffickToLys_Ab_Ag = - (p.ke*y(Ab_Ag_endo)) * (1 - p.frec); % Complex Trafficking from Endosome to Lysosome

% Ag: Production, Internalization, Recycling, Trafficking
Ag_Production = p.kAgProduction/p.Vmedia;
Int_Ag = - p.kiAg*y(Ag);%/y(Cells); % Antigen Internalization 
Rec_Ag = - (p.keAg*y(Ag_endo)) * p.frecAg;% / y(Cells); % Antigen Recycling 
TraffickToLys_Ag = - (p.keAg*y(Ag_endo)) * (1 - p.frecAg);% * y(Cells); % Complex Trafficking from Endosome to Lysosome

% ADC: Non-specific Uptake, Recycling, Trafficking
Int_ADC = - p.kiADC*(y(ADC)-y(ADC_endo)); % ADC Pinocytosis
Rec_ADC_FcRn = - p.krecADC*y(ADC_FcRn_endo);
TraffickToLys_ADC = - p.keADC*y(ADC_endo); % TraffickToLys_ADC = - p.ke*y(ADC_endo);

% Ab: Non-specific Uptake, Recycling, Trafficking
Int_Ab = - p.kiADC*(y(Ab)-y(Ab_endo)); % ADC Pinocytosis
Rec_Ab_FcRn = - p.krecADC*y(Ab_FcRn_endo); % FcRn Recycling
TraffickToLys_Ab = - p.keADC*y(Ab_endo); % TraffickToLys_Ab

% FcRn: Binding and Unbinding
BindUnbind_ADC_FcRn_endo = (- p.konFcRn*y(ADC_endo)*y(FcRn_endo) + p.koffFcRn*y(ADC_FcRn_endo)); % Intracellular Target-Antigen Binding and Unbinding
BindUnbind_Ab_FcRn_endo =  (- p.konFcRn*y(Ab_endo)*y(FcRn_endo) + p.koffFcRn*y(Ab_FcRn_endo)); % Intracellular Target-Antigen Binding and Unbinding

% Soluble Antigen
BindUnbind_ADC_Ag_s = - p.konAg_s*y(ADC)*y(Ag_s) + p.koffAg_s*y(ADC_Ag_s); % Target-Antigen Binding and Unbinding
BindUnbind_Ab_Ag_s =  - p.konAg_s*y(Ab)*y(Ag_s) + p.koffAg_s*y(Ab_Ag_s); % Target-Antigen Binding and Unbinding

% Intracellular Ag Binding and Unbinding
BindUnbind_ADC_Ag_endo = (- p.konAg*y(ADC_endo)*y(Ag_endo) + p.koffAg*y(ADC_Ag_endo)); % Intracellular Target-Antigen Binding and Unbinding
BindUnbind_Ab_Ag_endo = (- p.konAg*y(Ab_endo)*y(Ag_endo) + p.koffAg*y(Ab_Ag_endo)); % Intracellular Target-Antigen Binding and Unbinding

% Extracellular Degradation
ExtracellularDegradation_ADC = - p.kdeconjugation * y(ADC);

% Intracellular Degradation
IntracellularDegradation_ADC_Ag_lys = - y(ADC_Ag_lys);
IntracellularDegradation_ADC_lys = - y(ADC_lys);
IntracellularDegradation_Ag_lys = - y(Ag_lys);

IntracellularDegradation_Ab_lys = - y(Ab_lys);
IntracellularDegradation_Ab_Ag_lys = - y(Ab_Ag_lys);

% Warhead movement
Efflux_W = - p.keff * y(Win);
Influx_W = - p.kinf * y(Wex);% * InfluxVolumeCorrection;
InfluxVolumeCorrection = p.Vcytoplasm/p.Vmedia; % InfluxVolumeCorrection = (p.Vcytoplasm * y(Cells))/p.Vmedia;
TraffickToNuc_W = - p.kinN*y(Win); 
TraffickFromNuc_W = - p.koutN*y(Wn);
BindUnbind_Wn_DNA = - p.konDNA*y(Wn)*y(DNA) + p.koffDNA*y(Wn_DNA); % Represents DNA crosslink formation

% Volume Corrections
VolumeCorrection_Ex_to_Endo = p.Vmedia/p.Vendo;
VolumeCorrection_Endo_to_Ex = p.Vendo/p.Vmedia;
VolumeCorrection_Ex_to_Cyto = p.Vmedia/p.Vcytoplasm;
VolumeCorrection_Cyto_to_Ex = p.Vcytoplasm/p.Vmedia;
VolumeCorrection_Endo_to_Cyto = p.Vendo/p.Vcytoplasm;
VolumeCorrection_Cyto_to_Nuc = p.Vcytoplasm/p.Vnuc;
VolumeCorrection_Nuc_to_Cyto = p.Vnuc/p.Vcytoplasm;

%% Equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cell Surface
dydt(ADC) = p.q/p.Vmedia + y(Cells) * (BindUnbind_ADC_Ag + VolumeCorrection_Endo_to_Ex * Int_ADC - VolumeCorrection_Endo_to_Ex * Rec_ADC_FcRn)...
    + ExtracellularDegradation_ADC + BindUnbind_ADC_Ag_s;
dydt(Ag) =  Ag_Production + BindUnbind_ADC_Ag + BindUnbind_Ab_Ag + Int_Ag - VolumeCorrection_Endo_to_Ex * Rec_Ag;
dydt(ADC_Ag) =  - BindUnbind_ADC_Ag + Int_ADC_Ag - VolumeCorrection_Endo_to_Ex * Rec_ADC_Ag;

dydt(Ag_s) = BindUnbind_ADC_Ag_s + BindUnbind_Ab_Ag_s;
dydt(ADC_Ag_s) =  - BindUnbind_ADC_Ag_s;

dydt(Ab) = y(Cells) * (BindUnbind_Ab_Ag + VolumeCorrection_Endo_to_Ex * Int_Ab - VolumeCorrection_Endo_to_Ex * Rec_Ab_FcRn)...
    - ExtracellularDegradation_ADC + BindUnbind_Ab_Ag_s;
dydt(Ab_Ag) = - BindUnbind_Ab_Ag + Int_Ab_Ag - VolumeCorrection_Endo_to_Ex * Rec_Ab_Ag;
dydt(Ab_Ag_s) = - BindUnbind_Ab_Ag_s;


%% Endosome
dydt(ADC_endo) =     - Int_ADC + BindUnbind_ADC_Ag_endo + BindUnbind_ADC_FcRn_endo + TraffickToLys_ADC;
dydt(Ag_endo) =      - VolumeCorrection_Ex_to_Endo * Int_Ag + BindUnbind_ADC_Ag_endo + BindUnbind_Ab_Ag_endo + Rec_Ag + TraffickToLys_Ag;
dydt(ADC_Ag_endo) =  - VolumeCorrection_Ex_to_Endo * Int_ADC_Ag - BindUnbind_ADC_Ag_endo + Rec_ADC_Ag + TraffickToLys_ADC_Ag;

dydt(Ab_endo) =    - Int_Ab + BindUnbind_Ab_Ag_endo + BindUnbind_Ab_FcRn_endo + TraffickToLys_Ab;
dydt(Ab_Ag_endo) = - VolumeCorrection_Ex_to_Endo * Int_Ab_Ag - BindUnbind_Ab_Ag_endo + Rec_Ab_Ag + TraffickToLys_Ab_Ag;

dydt(FcRn_endo) =       BindUnbind_ADC_FcRn_endo - Rec_ADC_FcRn + BindUnbind_Ab_FcRn_endo - Rec_Ab_FcRn;
dydt(ADC_FcRn_endo) = - BindUnbind_ADC_FcRn_endo + Rec_ADC_FcRn;
dydt(Ab_FcRn_endo) =  - BindUnbind_Ab_FcRn_endo + Rec_Ab_FcRn;

%% Lysosome
dydt(ADC_lys) =    - TraffickToLys_ADC + IntracellularDegradation_ADC_lys;
dydt(Ag_lys) =     - TraffickToLys_Ag + IntracellularDegradation_Ag_lys;
dydt(ADC_Ag_lys) = - TraffickToLys_ADC_Ag + IntracellularDegradation_ADC_Ag_lys;

dydt(Ab_lys) =     - TraffickToLys_Ab + IntracellularDegradation_Ab_lys;
dydt(Ab_Ag_lys) =  - TraffickToLys_Ab_Ag + IntracellularDegradation_Ab_Ag_lys;

%% Clearance Compartment
dydt(ADC_cl) = - (p.Vendo * y(Cells)) * (IntracellularDegradation_ADC_lys + IntracellularDegradation_ADC_Ag_lys);
dydt(Ag_cl) = - (p.Vendo * y(Cells)) * (IntracellularDegradation_Ag_lys + IntracellularDegradation_ADC_Ag_lys + IntracellularDegradation_Ab_Ag_lys);
dydt(Ab_cl) = - (p.Vendo * y(Cells)) * (IntracellularDegradation_Ab_lys + IntracellularDegradation_Ab_Ag_lys);

%% Cytosol
dydt(Win) = - VolumeCorrection_Endo_to_Cyto  * p.DAR * ...
    (IntracellularDegradation_ADC_lys + IntracellularDegradation_ADC_Ag_lys)...
    + Efflux_W - VolumeCorrection_Ex_to_Cyto * Influx_W * InfluxVolumeCorrection + TraffickToNuc_W - VolumeCorrection_Nuc_to_Cyto * TraffickFromNuc_W;

% Intracellular Warhead Tracking
dydt(Win_lys) = - VolumeCorrection_Endo_to_Cyto * p.DAR * ...
    (IntracellularDegradation_ADC_lys + IntracellularDegradation_ADC_Ag_lys) - (p.keff + p.kinN) * y(Win_lys);
dydt(Win_deg) = VolumeCorrection_Ex_to_Cyto * p.kinf * y(Wex_deg) * InfluxVolumeCorrection - (p.keff + p.kinN) * y(Win_deg);
dydt(Win_influx) = VolumeCorrection_Ex_to_Cyto * p.kinf * y(Wex_efflux) * InfluxVolumeCorrection - (p.keff + p.kinN) * y(Win_influx);
dydt(Win_nuc) = - VolumeCorrection_Nuc_to_Cyto * TraffickFromNuc_W - (p.keff + p.kinN) * y(Win_nuc);

%% Extracellular Space
dydt(Wex) = - ExtracellularDegradation_ADC * p.DAR - VolumeCorrection_Cyto_to_Ex * Efflux_W * y(Cells) + Influx_W * InfluxVolumeCorrection * y(Cells);

% Extracellular Warhead Tracking
dydt(Wex_deg) = - ExtracellularDegradation_ADC * p.DAR - p.kinf * y(Wex_deg) * InfluxVolumeCorrection * y(Cells);
dydt(Wex_efflux) = VolumeCorrection_Cyto_to_Ex * p.keff * y(Win) * y(Cells) - p.kinf * y(Wex_efflux) * InfluxVolumeCorrection * y(Cells); %(y(Win_lys) + y(Win_deg) + y(Win_influx) + y(Win_nuc)) 

%% Nucleus
dydt(Wn) =       BindUnbind_Wn_DNA - VolumeCorrection_Cyto_to_Nuc * TraffickToNuc_W + TraffickFromNuc_W;
dydt(DNA) =      BindUnbind_Wn_DNA;
dydt(Wn_DNA) = - BindUnbind_Wn_DNA;

% Crosslink Formation
Crosslinks = 1 ./ (1 + (p.KA./y(Wn_DNA)).^p.n);

%% Cell Growth and Inhibition
dydt(Cells) = p.kg * y(Cells) - p.kkill * Crosslinks * y(Cells);


