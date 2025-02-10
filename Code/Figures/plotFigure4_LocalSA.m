doseType = "ADC"; % "ADC" or "Ab" or "PBD" or "IsotypeADC" or "None"
time = 24*4; % (hr) incubation time in hours
run(setup_file);

%% Univariate Sensitivity Analysis for different outputs
if plotFigure4A == "yes"

    fprintf("Generating Figure 4A...\n")

    % Pick a percent change for the parameter
    percentVaried = 10;
    SensVar = [Wex Win Wn_DNA Cells];
    SensVarNames = {'Warhead_e_x','Warhead_i_n','Warhead-DNA','Crosslinks','Cell Population'};
    allParamNames = fieldnames(p);
    excludedParams = {'q', 'kdeconjugation'}; % Parameters to exclude
    paramNames = setdiff(allParamNames, excludedParams,'stable'); % Remove excluded parameters from analysis
    numParams = numel(paramNames); 

    % Define variables
    T = cell(1,numParams);
    Y = cell(1,numParams);
    AUC = cell(numParams,numel(SensVar));
    Sens = cell(numParams,numel(SensVar));
    Crosslinks = cell(1,numParams);
    AUC_Crosslinks = cell(1,numParams);
    Sens_Crosslinks = cell(numParams,1);
    AUC0 = zeros(1,numel(SensVar));

    % Base Case
    y0 = InitCond;
    [T0,Y0] = ode23s(eqns_file,[0 time],y0,ode_options,p);
    for j = 1:numel(SensVar)
        AUC0(j) = trapz(T0,Y0(:,SensVar(j)));
    end
    Crosslinks0 = 1 ./ (1 + (p.KA./Y0(:,Wn_DNA)).^p.n); 
    AUC0_Crosslinks = trapz(T0,Crosslinks0);

    for i = 1:numParams % Loop through parameters in system

        pSensitivity = p;
        origParam = p.(paramNames{i});
        testParam = origParam * (1 + (percentVaried/100));
        pSensitivity.(paramNames{i}) = testParam;           

        [T{i},Y{i}] = ode23s(eqns_file,[0 time],y0,ode_options,pSensitivity);

        for j = 1:numel(SensVar)
            AUC{i,j} = trapz(T{i},Y{i}(:,SensVar(j))); % AUC of Selected Var
            Sens{i,j} = ((AUC{i,j}-AUC0(j))/AUC0(j))/((testParam-origParam)/origParam);
        end

        % For Crosslinking
        Crosslinks{i} = 1 ./ (1 + (pSensitivity.KA./Y{i}(:,Wn_DNA)).^pSensitivity.n); 
        AUC_Crosslinks{i} = trapz(T{i},Crosslinks{i}); % AUC of Selected Var
        Sens_Crosslinks{i,1} = ((AUC_Crosslinks{i}-AUC0_Crosslinks)/AUC0_Crosslinks)/((testParam-origParam)/origParam);

        disp(['Finished loop ',num2str(i),' out of ',num2str(numParams)])            

    end

    SensMatrix = cell2mat(Sens);
    SensMatrix_Crosslinks = cell2mat(Sens_Crosslinks);

    SensMatrix_Final = horzcat(SensMatrix(:,1:3),SensMatrix_Crosslinks,SensMatrix(:,4));

    fig = figure();
    clim(gca,[min(SensMatrix_Final,[],"all") max(SensMatrix_Final,[],"all")]);
    HM = heatmap(SensVarNames,paramNames,SensMatrix_Final,...
        'FontSize',20,'Colormap',bluewhitered,'CellLabelColor','none');
    set(gca,'FontSize',20)
    set(gcf,'color','w','position',[500 500 750 1000])
    hs = struct(HM);
    ylabel(hs.Colorbar, "Sensitivity",'FontSize',20);
    drawnow;
    fprintf("Done!\n")

end

% Univariate local sensitivity analysis at different time points
if plotFigure4B == "yes" 

    fprintf("Generating Figure 4B...\n")

    % Pick a percent change for the parameter
    percentVaried = 10;
    
    % Time Cases
    titles = {'1 hour',...
              '12 hours',...
              '24 hours'...                 
              '48 hours',...
              '72 hours',...
              '96 hours','168 hours'};
    time_Case = numel(titles);
    time_values = [1 12 24 48 72 96 168];
    allParamNames = fieldnames(p);
    excludedParams = {'q', 'kdeconjugation'}; % Parameters to exclude
    paramNames = setdiff(allParamNames, excludedParams,'stable'); % Remove excluded parameters from analysis
    numParams = numel(paramNames);
    T = cell(time_Case,numParams);
    Y = cell(time_Case,numParams);
    AUC = cell(time_Case,numParams);
    Sens = cell(time_Case,numParams);
    SensVar = Cells;

    for i = 1:time_Case % Loop through cases
       
        % Base Case
        y0 = InitCond;
        [T0,Y0] = ode23s(eqns_file,[0 time_values(i)],y0,ode_options,p);
        AUC0 = trapz(T0,Y0(:,SensVar));

        for j = 1:numParams % Loop through parameters in system
        
            pSensitivity = p;
            origParam = p.(paramNames{j});
            testParam = origParam * (1 + (percentVaried/100));
            pSensitivity.(paramNames{j}) = testParam;           

            [T{i,j},Y{i,j}] = ode23s(eqns_file,[0 time_values(i)],y0,ode_options,pSensitivity);
            AUC{i,j} = trapz(T{i,j},Y{i,j}(:,SensVar));             

            % Calculate Sensitivity
            Sens{i,j} = ((AUC{i,j}-AUC0)/AUC0)/((testParam-origParam)/origParam);                
        end

        disp(['Finished loop ',num2str(i),' out of ',num2str(time_Case)])
        
    end
    
    SensMatrix = cell2mat(Sens)';

    fig = figure();
    clim(gca,[min(SensMatrix,[],"all") max(SensMatrix,[],"all")]);
    HM = heatmap(titles,paramNames,SensMatrix,...
        'FontSize',20,'Colormap',bluewhitered,'CellLabelColor','none');
    set(gca,'FontSize',20)
    set(gcf,'color','w','position',[500 500 750 1000])
    hs = struct(HM);
    ylabel(hs.Colorbar, "Sensitivity",'FontSize',20);
    drawnow;
    fprintf("Done!\n")

end

% Univariate local sensitivity analysis at different receptor expression
if plotFigure4C == "yes" 

    fprintf("Generating Figure 4C...\n")

    % Pick a percent change for the parameter
    percentVaried = 10;
    
    % Receptor Cases
    titles = {'10^1','10^2',...
              '10^3',...
              '10^4'...                 
              '10^5',...
              '10^6','10^7'};
    receptorCase = numel(titles);
    
    allParamNames = fieldnames(p);
    excludedParams = {'q', 'kdeconjugation'}; % Parameters to exclude
    paramNames = setdiff(allParamNames, excludedParams,'stable'); % Remove excluded parameters from analysis
    numParams = numel(paramNames);
    T = cell(receptorCase,numParams);
    Y = cell(receptorCase,numParams);
    AUC = cell(receptorCase,numParams);
    Sens = cell(receptorCase,numParams);
    InitConds = repmat(InitCond,1,receptorCase);
    InitConds(Ag,:) = logspace(1,receptorCase,receptorCase) / NA * 10^9 / p.Vmedia;
    
    % Variable to perform sensitivity on
    SensVar = Cells;

    for i = 1:receptorCase % Loop through cases
       
        % Base Case
        y0 = InitConds(:,i);
        [T0,Y0] = ode23s(eqns_file,[0 time],y0,ode_options,p);
        AUC0 = trapz(T0,Y0(:,SensVar)); % AUC of Cells

        for j = 1:numParams % Loop through parameters in system
        
            pSensitivity = p;
            origParam = p.(paramNames{j});
            testParam = origParam * (1 + (percentVaried/100));
            pSensitivity.(paramNames{j}) = testParam;           
            [T{i,j},Y{i,j}] = ode23s(eqns_file,[0 time],y0,ode_options,pSensitivity);
            AUC{i,j} = trapz(T{i,j},Y{i,j}(:,SensVar)); % AUC of Cells
            Sens{i,j} = ((AUC{i,j}-AUC0)/AUC0)/((testParam-origParam)/origParam);
            
        end

        disp(['Finished loop ',num2str(i),' out of ',num2str(receptorCase)])
        
    end
    
    SensMatrix = cell2mat(Sens)';

    fig = figure();
    clim(gca,[min(SensMatrix,[],"all") max(SensMatrix,[],"all")]);
    HM = heatmap(titles,paramNames,SensMatrix,...
        'FontSize',20,'Colormap',bluewhitered,'CellLabelColor','none');
    set(gca,'FontSize',20)
    set(gcf,'color','w','position',[500 500 750 1000])
    hs = struct(HM);
    ylabel(hs.Colorbar, "Sensitivity",'FontSize',20);
    xlabel('\bfReceptors Per Cell')
    drawnow;
    fprintf("Done!\n")

end

if plotFigure4D == "yes" 

    fprintf("Generating Figure 4D...\n")
    
    % Pick a percent change for the parameter
    percentVaried = 10;
    
    % Dose Cases
    titles = {'10^-^5',...
              '10^-^4',...
              '10^-^3'...                 
              '10^-^2',...
              '10^-^1',...
              '10^0',...
              '10^1'};
    doseCase = numel(titles);
    dose_values = [1e-5 1e-4 1e-3 1e-2 1e-1 1e0 1e1]; %logspace(-5,1,doseCase);
    
    allParamNames = fieldnames(p);
    excludedParams = {'q', 'kdeconjugation'}; % Parameters to exclude
    paramNames = setdiff(allParamNames, excludedParams,'stable'); % Remove excluded parameters from analysis
    numParams = numel(paramNames);
    T = cell(doseCase,numParams);
    Y = cell(doseCase,numParams);
    AUC = cell(doseCase,numParams);
    Sens = cell(doseCase,numParams);
    
    % Variable to perform sensitivity on
    SensVar = Cells;

    for i = 1:doseCase % Loop through cases
       
        % Base Case
        y0 = InitCond;
        y0(ADC) = dose_values(i);
        [T0,Y0] = ode23s(eqns_file,[0 time],y0,ode_options,p);
        AUC0 = trapz(T0,Y0(:,SensVar)); % AUC of Cells

        for j = 1:numParams % Loop through parameters in system
        
            pSensitivity = p;
            origParam = p.(paramNames{j});
            testParam = origParam * (1 + (percentVaried/100));
            pSensitivity.(paramNames{j}) = testParam;           
            [T{i,j},Y{i,j}] = ode23s(eqns_file,[0 time],y0,ode_options,pSensitivity);
            AUC{i,j} = trapz(T{i,j},Y{i,j}(:,SensVar)); % AUC of Cells
            Sens{i,j} = ((AUC{i,j}-AUC0)/AUC0)/((testParam-origParam)/origParam);
            
        end

        disp(['Finished loop ',num2str(i),' out of ',num2str(doseCase)])
        
    end
    
    SensMatrix = cell2mat(Sens)';

    fig = figure();
    clim(gca,[min(SensMatrix,[],"all") max(SensMatrix,[],"all")]);
    HM = heatmap(titles,paramNames,SensMatrix,...
        'FontSize',20,'Colormap',bluewhitered,'CellLabelColor','none');
    set(gca,'FontSize',20)
    set(gcf,'color','w','position',[500 500 750 1000])
    hs = struct(HM);
    ylabel(hs.Colorbar, "Sensitivity",'FontSize',20);
    xlabel('\bfADC Dose (nM)')
    drawnow;
    fprintf("Done!\n")

end