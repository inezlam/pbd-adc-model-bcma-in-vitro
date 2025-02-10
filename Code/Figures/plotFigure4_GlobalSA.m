doseType = "ADC"; % "ADC" or "Ab" or "PBD" or "IsotypeADC" or "None"
time = 24*4; % (hr) incubation time in hours
run(setup_file);

%% Global Sensitivity Analysis %%%%%%%%%%%%%%%%%%%
if plotFigure4E == "yes"

    fprintf("Generating Figure 4E...\n")

    parameterMultiplier = logspace(-2,2,9);
    allParamNames = fieldnames(p);
    excludedParams = {'Vmedia', 'Vendo', 'Vnuc', 'Vcytoplasm', 'q', 'DAR', 'frec', 'frecAg', 'n'}; % Parameters to exclude
    paramNames = setdiff(allParamNames, excludedParams,'stable'); % Remove excluded parameters from analysis
    numParams = numel(paramNames); 
    globalSens = cell(numParams,numel(parameterMultiplier));
    
    
    ADC_doses = [1.85999E-05 5.58656E-05 0.000167531 0.000502644 0.001507935 0.004523797 0.013571392 0.040714182 0.122142546 0.366427639 1.099283048 3.297848484]; % nM - initial concentrations of ADC        
    x_exp = ADC_doses;
    y_exp = [92.33 103.55 97.36 106.42 98.37 102.4 70.33 50.05 29.77 16.68 6.62 2.44]; % ONC2228-0004 Fig 13.1-1 In vitro cytotox_cell panel
    
    for i = 1:numParams % Loop through cases
    
        pTest = p;
       
        for j = 1:numel(parameterMultiplier) % Loop through parameters in system
            
            % Change parameter to test with multiplier
            pTest.(paramNames{i}) = p.(paramNames{i})*parameterMultiplier(j);       
    
            % Calculate error for global sensitivity
            globalSens{i,j} = rssq(GS_costfxn(eqns_file, x_exp, y_exp, pTest, time, ADC_doses, ADC, InitCond));
                  
            % disp(['   Finished loop ',num2str(j),' out of ',num2str(numel(parameterMultiplier))])            
            
        end
    
        disp(['Finished loop ',num2str(i),' out of ',num2str(numParams)])
        
    end
    
    GlobalSensMatrix = cell2mat(globalSens);

    multiplierNames = ["-2"," ","-1"," ","0"," ","1"," ","2"];

    fig = figure();
    GS_HM = heatmap(linspace(-2,2,9),paramNames,GlobalSensMatrix,...
        'FontSize',20,'Colormap',flip(pink),'CellLabelColor','none');
    GS_HM.XDisplayLabels = multiplierNames;
    xlabel('\bfChange in Orders of Magnitude (log_1_0)')
    set(gcf,'color','w','position',[500 500 700 750])
    hs = struct(GS_HM);
    ylabel(hs.Colorbar, "Cost",'FontSize',20);
    drawnow;
    fprintf("Done!\n")

end