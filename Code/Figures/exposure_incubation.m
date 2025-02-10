function [T,Y] = exposure_incubation(fxn, p, ExposureTime, IncubationTime, InitCond, DrugSpecies)
%{ 
Drug Exposure followed by Incubation Simulation
T = time
Y = concentration of molecular species
fxn = equations file used
p = parameters
InitCond = initial conditions
ExposureTime = time exposed to drug
IncubationTime = time inbuated following exposure
%}

%     fxn=@BCMA_eqns_crosslink;
%     ExposureTime=2;
%     IncubationTime=24;
%     InitCond = InitCond_PBD;

%     fprintf('Simulating drug exposure and incubation...')
    ode_options = odeset('MaxStep',5e-1, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);

    %% Simulate Exposure Period
    [T1,Y1] = ode23s(fxn,[0 ExposureTime],InitCond,ode_options,p);        

    %% Simulate Incubation Period
    IncubationInitCond = Y1(end,:); % Ending of exposure simulation becomes IC for incubation simulation
    IncubationInitCond(DrugSpecies) = 0; % Existing extracelluar drug is washed out
    [T2,Y2] = ode23s(fxn,[0 IncubationTime],IncubationInitCond,ode_options,p);
    
    %% Combine and return
    T = vertcat(T1(1:end-1),T2+ExposureTime);
    Y = vertcat(Y1(1:end-1,:),Y2);

%     fprintf(' Done!\n')
end 