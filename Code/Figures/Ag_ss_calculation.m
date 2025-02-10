%{
Function to calculate Antigen levels at steady state

Inputs:
- fxn = equations file to use
- p = structure containing parameters
- time = time used for steady state simulation
- AgSpecies = numbering for Antigen species
- InitCond = initial conditions for simulation
- AgProd = antigen production rate used
- options = ODE solver options

Outputs: Ag_ss (Concentration of Ag at steady state)
%}

function [T, Ag_ss, Ag_endo_ss, Ag_lys_ss] = Ag_ss_calculation2(fxn,ode_options,p,time,AgSpecies,AgEndoSpecies,AgLysSpecies,InitCond,AgProd)
  
    % Set up steady state simulation
    p_temp = p;
    p_temp.kg = 0;
    p_temp.kAgProduction = AgProd;
%     [~,Y0] = ode23s(fxn,[0 time],InitCond,ode_options,p_temp);   
    [T,Y] = ode45(fxn,[0 time],InitCond,ode_options,p_temp);   

    % Calculate Ag level at steady state
    Ag_ss = Y(:,AgSpecies);
    Ag_endo_ss = Y(:,AgEndoSpecies);
    Ag_lys_ss = Y(:,AgLysSpecies);

    
end