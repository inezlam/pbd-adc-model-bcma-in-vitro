function [X,crosslink_percentage] = dose_response_crosslink_incubation(fxn, p, ExposureTime, IncubationTime, doses, DrugSpecies, CrosslinkSpecies, InitCond)
%{
  X = doses
  crosslink_percentage = output
  fxn = equations file
  p = parameters
  ExposureTime = time of exposure to drug
  IncubationTime = time of incubation without drug
  doses = doses given of drug
  Species = molecular species of interest
  InitCond = initial conditions
%}

%     fprintf('Calculating dose response...')
%     options = odeset('MaxStep',5e-1, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
    
    X = doses;
    crosslink_percentage = nan(1,length(doses)); % Crosslink percentage after incubation time
    
    % Calculate population of untreated control cells
%     InitCond(Species) = 0;
%     [~,Y0] = ode23s(fxn,[0 time],InitCond,options,p);        
%     untreated_cell_population = Y0(end,length(InitCond));    
    
    % Calculated population of treated cells
%     figure
%     hold on
    for dose = 1:length(doses)
        InitCond(DrugSpecies) = doses(dose); % nM   
        [T,Y] = exposure_incubation(fxn, p, ExposureTime, IncubationTime, InitCond, DrugSpecies);
%         [T,Y] = ode23s(fxn,[0 time],InitCond,options,p);   

%         plot(T,Y(:,end))
        Crosslinks = 1 ./ (1 + (p.KA./Y(:,CrosslinkSpecies)).^p.n);
        
%         [~,closestIndex] = min(abs(T-IncubationTime));
%         Crosslink_Norm = 0.65/Crosslinks(closestIndex);
        
%         plot(T,Crosslinks)
        crosslink_percentage(dose) = Crosslinks(end);
    end
    
%     crosslink_percentage = 100*(crosslink_percentage./Y(end-2,1));
%     cell_survival = cell_survival/cell_survival(1,1) * 100;
%     cell_survival = cell_population*100/InitCond(end);
%     Y_end = Y_end*100/max(Y_end);
%     fprintf(' Done!\n')
end 