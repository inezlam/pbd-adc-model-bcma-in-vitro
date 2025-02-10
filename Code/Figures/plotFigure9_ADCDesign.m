doseType = "ADC"; % "ADC" or "Ab" or "PBD" or "IsotypeADC" or "None"
time = 24*4; % (hr) incubation time in hours
run(setup_file);

if plotFigure9 == "yes"

    fprintf("Generating Figure 9...\n")

    %% Trivariate Analysis - Heatmaps
    DAR_values = [1 2 3 4 5 6 7 8];
    kkill_values = [1e-3 1e-2 1e-1 1e0 1e1 1e2]; 
    R_values = logspace(-1,3,5);

    SensVar1 = Cells;
    SensVar2 = Win_influx;
    SensVar3 = Win_influx;

    for i = 1:numel(DAR_values)

        pTest = p;
        pTest.DAR = DAR_values(i);

        for j = 1:numel(kkill_values) % Loop through parameters in system
        
            pTest.kkill = kkill_values(j);

            for k = 1:numel(R_values)

                pTest.keff = pTest.kinf/(1 + R_values(k));
                [T{i,j,k},Y{i,j,k}] = ode23s(eqns_file,[0 time],InitCond,ode_options,pTest);
                [T0{i,j,k},Y0{i,j,k}] = ode23s(eqns_file,[0 time],InitCond_noADC,ode_options,pTest);
                CellSurvival{i,j,k} = Y{i,j,k}(end,Cells)/Y0{i,j,k}(end,Cells)*100;
                AUC1{i,j,k} = trapz(T{i,j,k},Y{i,j,k}(:,SensVar1)); % AUC 
                AUC2{i,j,k} = trapz(T{i,j,k},Y{i,j,k}(:,SensVar2)); % AUC 
                AUC3{i,j,k} = trapz(T{i,j,k},Y{i,j,k}(:,SensVar3)) .* pTest.kkill; % AUC 
           end
           
        end

        disp(['Finished loop ',num2str(i),' out of ',num2str(numel(DAR_values))])   
        
    end

    CellSurvivalMatrix = cell2mat(CellSurvival);
    AUCMatrix1 = cell2mat(AUC1);
    AUCMatrix2 = cell2mat(AUC2);
    AUCMatrix3 = cell2mat(AUC3);

    %% Plots
    fig = figure();
    t = tiledlayout(3,numel(R_values),'TileSpacing','Compact','TileIndexing','columnmajor');
    set(gcf,'color','w','position',[100 100 1700 900])

    for r = 1:numel(R_values)
        % Cell Population - AUC
        nexttile
        HM2 = heatmap(DAR_values,kkill_values,AUCMatrix1(:,:,r)',...
            'FontSize',20,'Colormap',bluewhitered,'CellLabelColor','none');
        HM2.Title = ['\bf R = ' num2str(R_values(r))];
        HM2.XDisplayLabels = nan(size(DAR_values));
        if r == 1
        elseif r > 1
            HM2.YDisplayLabels = nan(size(kkill_values));
        end
        if r < numel(R_values)
            HM2.ColorbarVisible = 'off';
        end
        set(gca,'FontSize',20)
        clim(gca,[min(AUCMatrix1,[],"all") max(AUCMatrix1,[],"all")]);
        hs2 = struct(HM2);
        ylabel(hs2.Colorbar, ["Cell Population (AUC)"],'FontSize',20,'FontWeight','bold');

        % Influxed Warhead - AUC
        nexttile
        HM3 = heatmap(DAR_values,kkill_values,AUCMatrix2(:,:,r)',...
            'FontSize',20,'Colormap',bluewhitered,'CellLabelColor','none');
        HM3.XDisplayLabels = nan(size(DAR_values));
        if r == 1
            HM3.YLabel = '\bf k_{kkill}'; 
        elseif r > 1
            HM3.YDisplayLabels = nan(size(kkill_values));
        end
        if r < numel(R_values)
            HM3.ColorbarVisible = 'off';
        end
        set(gca,'FontSize',20)
        clim(gca,[min(AUCMatrix2,[],"all") max(AUCMatrix2,[],"all")]);
        hs3 = struct(HM3);
        ylabel(hs3.Colorbar, ["Influxed Warhead (AUC)"],'FontSize',20,'FontWeight','bold');

        % Influxed Warhead * kkill
        nexttile
        HM4 = heatmap(DAR_values,kkill_values,AUCMatrix3(:,:,r)',...
            'FontSize',20,'Colormap',bluewhitered,'CellLabelColor','none');
        if r == 1
        elseif r > 1
            HM4.YDisplayLabels = nan(size(kkill_values));
        end
        if r == 3
            HM4.XLabel = '\bf DAR'; 
        end
        if r < numel(R_values)
            HM4.ColorbarVisible = 'off';
        end
        set(gca,'FontSize',20)
        clim(gca,[min(AUCMatrix3,[],"all") max(AUCMatrix3,[],"all")]);
        hs4 = struct(HM4);
        ylabel(hs4.Colorbar, ["Influxed Warhead * kkill"],'FontSize',20,'FontWeight','bold');

    end

    drawnow;
    fprintf("Done!\n")
end