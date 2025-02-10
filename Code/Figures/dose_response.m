function [X,cell_survival] = dose_response(fxn, p, time, doses, DoseSpecies, InitCond)

%     fprintf('Calculating dose response...')
    ode_options = odeset('MaxStep',5e-1, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);
    
    X = doses;
    cell_population = nan(1,length(doses)); % Cell Survival after incubation time
    
    % Calculate population of untreated control cells
    InitCond(DoseSpecies) = 0;
    [~,Y0] = ode23s(fxn,[0 time],InitCond,ode_options,p);        
    untreated_cell_population = Y0(end,length(InitCond));    
    
    % Calculated population of treated cells
%     figure
%     hold on
    for dose = 1:length(doses)
        InitCond(DoseSpecies) = doses(dose); % nM   
        [T,Y] = ode23s(fxn,[0 time],InitCond,ode_options,p);        
%         plot(T,Y(:,end))
        cell_population(dose) = Y(end,length(InitCond));
    end
    
    cell_survival = 100*(cell_population/untreated_cell_population);
%     cell_survival = cell_survival/cell_survival(1,1) * 100;
%     cell_survival = cell_population*100/InitCond(end);
%     Y_end = Y_end*100/max(Y_end);
%     fprintf(' Done!\n')
end            