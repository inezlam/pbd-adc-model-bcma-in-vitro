doseType = "ADC"; % "ADC" or "Ab" or "PBD" or "IsotypeADC" or "None"
time = 24*7; % (hr) incubation time in hours
run(setup_file);

if plotFigure8A == "yes"
    
    fprintf("Generating Figure 8A...\n")

    ExtracellularVolumes = logspace(-8,-4,5)';

    % For a range of volumes from Vmedia = Vcytoplasm to Vmedia >> Vcytoplasm
    for i = 1:numel(ExtracellularVolumes)
    
        pTest = p;
        pTest.Vmedia = ExtracellularVolumes(i);
        InitCond(Ag) = Ag0/pTest.Vmedia; % Initial extracellular membrane-bound Ag 
        InitCond_noADC(Ag) = Ag0/pTest.Vmedia;

        % Simulate treated cells
        [T{i},Y{i}] = ode23s(eqns_file,[0 time],InitCond,ode_options,pTest);

        % Simulate untreated cells (control)
        [T_control{i},Y_control{i}] = ode23s(eqns_file,T{i},InitCond_noADC,ode_options,pTest);        
  
        disp(['Finished loop ',num2str(i),' out of ',num2str(numel(ExtracellularVolumes))])
    
    end

    figRows = 3;
    figCols = numel(ExtracellularVolumes);
    fig8a = figure;
    tiledlayout(figRows,figCols,'TileSpacing','tight')
    set(gcf,'color','w','position',[200 200 1400 800])
    ylim_UB_Wex = zeros(numel(ExtracellularVolumes),1);
    ylim_UB_Win = zeros(numel(ExtracellularVolumes),1);

    for i = 1:numel(ExtracellularVolumes)
        ylim_UB_Wex(i,1) = max(Y{i}(:,Wex));
        ylim_UB_Win(i,1) = max(Y{i}(:,Win));
    end
    
    % Extracellular Warhead
    for i = 1:numel(ExtracellularVolumes)
        nexttile 
        area(T{i},Y{i}(:,Wex),'linewidth',3,'FaceColor','k','FaceAlpha',.1)
        hold on
        area(T{i},Y{i}(:,[Wex_deg, Wex_efflux]),'linewidth',3,'FaceAlpha',0.9)
        box on; grid on;
        xlim([0 time])
        ylim([0 max(ylim_UB_Wex)*(11/10)])
        set(gca,'XTickLabel',[],'FontSize',20)
        title(['10^-^' num2str(abs(log10(ExtracellularVolumes(i)))) ' L'],'fontsize',22)
        if i == 1
            ylabel({'Extracellular','Warhead (nM)'},'FontWeight','bold','fontsize',22)
        elseif i == numel(ExtracellularVolumes)
            h = legend({'',['Extracellular' newline 'Deconjugation'],'Efflux'},'Location','eastoutside');
            title(h,'Source')   
        end
        if i > 1
            set(gca,'YTickLabel',[])
        end
    end
    
    % Intracellular Warhead
    for i = 1:numel(ExtracellularVolumes)
        nexttile 
        area(T{i},Y{i}(:,Win),'linewidth',3,'FaceColor','k','FaceAlpha',.1)
        hold on
        area(T{i},Y{i}(:,[Win_deg, Win_influx, Win_lys, Win_nuc]),'linewidth',3,'FaceAlpha',1)
        box on; grid on;
        xlim([0 time])
        ylim([0 max(ylim_UB_Win)*(11/10)])
        set(gca,'XTickLabel',[],'FontSize',20)
        if i == 1
            ylabel({'Intracellular','Warhead (nM)'},'FontWeight','bold','fontsize',22)
        elseif i == numel(ExtracellularVolumes)
            h = legend({'',['Extracellular' newline 'Deconjugation'],'Influx','Lysosome','Nucleus'},'Location','eastoutside');
            title(h,'Source')
        end
        if i > 1
            set(gca,'YTickLabel',[])
        end
    end

    % Cells
    for i = 1:numel(ExtracellularVolumes)
        nexttile 
        plot(T{i},Y{i}(:,Cells)./Y_control{i}(:,Cells)*100,'linewidth',3,'Color','k')
        box on; grid on;
        xlim([0 time])
        ylim([0 100])
        set(gca,'FontSize',20)
        if i == 1
            ylabel({'Cell Survival','(%)'},'FontWeight','bold','fontsize',22)
        elseif i == 3
            xlabel('Time (hours)','fontsize',24,'FontWeight','bold')
        elseif i == numel(ExtracellularVolumes)
        end
        if i > 1
            set(gca,'YTickLabel',[])
        end
    end

    drawnow;
    fprintf("Done!\n")


    % %% Stacked Bar Plot
    % Xex = categorical([varNames(Wex),strcat(varNames(Wex_deg),{' + '},varNames(Wex_efflux))]);
    % Xex = reordercats(Xex,[varNames(Wex),strcat(varNames(Wex_deg),{' + '},varNames(Wex_efflux))]);
    % Xin = categorical([varNames(Win),strcat(varNames(Win_deg),{' + '},varNames(Win_influx),{' + '},varNames(Win_lys),{' + '},varNames(Win_nuc))]);
    % Xin = reordercats(Xin,[varNames(Win),strcat(varNames(Win_deg),{' + '},varNames(Win_influx),{' + '},varNames(Win_lys),{' + '},varNames(Win_nuc))]);
    % figRows = 2;
    % figCols = numel(ExtracellularVolumes);
    % figure; 
    % tiledlayout(figRows,figCols,'TileIndexing','columnmajor','TileSpacing','tight')
    % set(gcf,'color','w','position',[200 200 1600 950])
    % ylim_UB_WexAUC = zeros(numel(ExtracellularVolumes),1);
    % ylim_UB_WinAUC = zeros(numel(ExtracellularVolumes),1);
    % 
    % for i = 1:numel(ExtracellularVolumes)
    %     AUC1{i} = trapz(T{i},Y{i});
    %     ylim_UB_WexAUC(i,1) = max(AUC1{i}(Wex));
    %     ylim_UB_WinAUC(i,1) = max(AUC1{i}(Win));
    % end
    % 
    % for i = 1:numel(ExtracellularVolumes)
    % 
    %     nexttile 
    %     bar(AUC1{i}(Wex),'FaceColor','k','FaceAlpha',.1)
    %     hold on
    %     bar([ 1 ;nan],[AUC1{i}(Wex_deg) AUC1{i}(Wex_efflux); nan(1,2)],'stacked','FaceColor','flat');      
    %     hold off
    %     if i == 1
    %         ylabel('Extracellular Warhead (AUC)','fontsize',20,'FontWeight','bold')
    %     end
    %     if i > 1
    %         set(gca,'YTickLabel',[])
    %     end
    %     if i == numel(ExtracellularVolumes)
    %         h = legend({'',['Extracellular' newline 'Deconjugation'],'Efflux'},'Location','eastoutside');
    %         title(h,'Source')
    %     end
    %     ylim([0 0.22])
    %     set(gca,'FontSize',20,'XTickLabel',[])
    %     title(['10^-^' num2str(abs(log10(ExtracellularVolumes(i)))) ' L'],'fontsize',20,'FontWeight','bold')
    %     grid on; box on;
    % 
    %     nexttile 
    %     bar(AUC1{i}(Win),'FaceColor','k','FaceAlpha',.1)
    %     hold on
    %     bar([ 1 ;nan],[AUC1{i}(Win_deg) AUC1{i}(Win_influx) AUC1{i}(Win_lys) AUC1{i}(Win_nuc); nan(1,4)],'stacked','FaceColor','flat');      
    %     hold off
    %     if i == 1
    %         ylabel('Intracellular Warhead (AUC)','fontsize',20,'FontWeight','bold')
    %     end
    %     if i > 1
    %         set(gca,'YTickLabel',[])
    %     end
    %     if i == numel(ExtracellularVolumes)
    %         h = legend({'',['Extracellular' newline 'Deconjugation'],'Influx','Lysosome','Nucleus'},'Location','eastoutside');
    %         title(h,'Source')
    %     end
    %     ylim([0 20])
    %     set(gca,'FontSize',20,'XTickLabel',[])
    %     grid on; box on;
    % end
    % drawnow;
    % fprintf("Done!\n")
    
end

if plotFigure8B == "yes"

    fprintf("Generating Figure 8B...\n")

    p.kdeconjugation = kdeconjugation_MouseSerum;
    ExtracellularVolumes = logspace(-8,-4,5)';

    % For a range of volumes from Vmedia = Vcytoplasm to Vmedia >> Vcytoplasm
    for i = 1:numel(ExtracellularVolumes)
    
        pTest = p;
        pTest.Vmedia = ExtracellularVolumes(i);
        InitCond(Ag) = Ag0/pTest.Vmedia; % Initial extracellular membrane-bound Ag 
        InitCond_noADC(Ag) = Ag0/pTest.Vmedia;

        % Simulate treated cells
        [T{i},Y{i}] = ode23s(eqns_file,[0 time],InitCond,ode_options,pTest);

        % Simulate untreated cells (control)
        [T_control{i},Y_control{i}] = ode23s(eqns_file,T{i},InitCond_noADC,ode_options,pTest);        
   
        disp(['Finished loop ',num2str(i),' out of ',num2str(numel(ExtracellularVolumes))])
    
    end

    figRows = 3;
    figCols = numel(ExtracellularVolumes);
    
    fig8b = figure;
    tiledlayout(figRows,figCols,'TileSpacing','tight')
    set(gcf,'color','w','position',[200 200 1400 800])
    ylim_UB_Wex = zeros(numel(ExtracellularVolumes),1);
    ylim_UB_Win = zeros(numel(ExtracellularVolumes),1);

    for i = 1:numel(ExtracellularVolumes)
        ylim_UB_Wex(i,1) = max(Y{i}(:,Wex));
        ylim_UB_Win(i,1) = max(Y{i}(:,Win));
    end

    % Extracellular Warhead
    for i = 1:numel(ExtracellularVolumes)
        nexttile 
        area(T{i},Y{i}(:,Wex),'linewidth',3,'FaceColor','k','FaceAlpha',.1)
        hold on
        area(T{i},Y{i}(:,[Wex_deg, Wex_efflux]),'linewidth',3,'FaceAlpha',0.9)
        box on; grid on;
        xlim([0 time])
        ylim([0 max(ylim_UB_Wex)*(11/10)])
        set(gca,'XTickLabel',[],'FontSize',20)
        title(['10^-^' num2str(abs(log10(ExtracellularVolumes(i)))) ' L'],'fontsize',22)
        if i == 1
            ylabel({'Extracellular','Warhead (nM)'},'FontWeight','bold','fontsize',22)
        elseif i == numel(ExtracellularVolumes)
            h = legend({'',['Extracellular' newline 'Deconjugation'],'Efflux'},'Location','eastoutside');
            title(h,'Source')    
        end
        if i > 1
            set(gca,'YTickLabel',[])
        end
    end
    
    % Intracellular Warhead
    for i = 1:numel(ExtracellularVolumes)
        nexttile 
        area(T{i},Y{i}(:,Win),'linewidth',3,'FaceColor','k','FaceAlpha',.1)
        hold on
        area(T{i},Y{i}(:,[Win_deg, Win_influx, Win_lys, Win_nuc]),'linewidth',3,'FaceAlpha',1)
        box on; grid on;
        xlim([0 time])
        ylim([0 max(ylim_UB_Win)*(11/10)])
        set(gca,'XTickLabel',[],'FontSize',20)
        if i == 1
            ylabel({'Intracellular','Warhead (nM)'},'FontWeight','bold','fontsize',22)
        elseif i == numel(ExtracellularVolumes)
            h = legend({'',['Extracellular' newline 'Deconjugation'],'Influx','Lysosome','Nucleus'},'Location','eastoutside');
            title(h,'Source')
        end
        if i > 1
            set(gca,'YTickLabel',[])
        end
    end

    % Cells
    for i = 1:numel(ExtracellularVolumes)
        nexttile 
        plot(T{i},Y{i}(:,Cells)./Y_control{i}(:,Cells)*100,'linewidth',3,'Color','k')
        box on; grid on;
        xlim([0 time])
        ylim([0 100])
        set(gca,'FontSize',20)
        if i == 1
            ylabel({'Cell Survival','(%)'},'FontWeight','bold','fontsize',22)
        elseif i == 3
            xlabel('Time (hours)','fontsize',24,'FontWeight','bold')
        elseif i == numel(ExtracellularVolumes)
        end
        if i > 1
            set(gca,'YTickLabel',[])
        end
    end

    drawnow;
    fprintf("Done!\n")

    % % Stacked Bar Plot
    % Xex = categorical([varNames(Wex),strcat(varNames(Wex_deg),{' + '},varNames(Wex_efflux))]);
    % Xex = reordercats(Xex,[varNames(Wex),strcat(varNames(Wex_deg),{' + '},varNames(Wex_efflux))]);
    % Xin = categorical([varNames(Win),strcat(varNames(Win_deg),{' + '},varNames(Win_influx),{' + '},varNames(Win_lys),{' + '},varNames(Win_nuc))]);
    % Xin = reordercats(Xin,[varNames(Win),strcat(varNames(Win_deg),{' + '},varNames(Win_influx),{' + '},varNames(Win_lys),{' + '},varNames(Win_nuc))]);
    % 
    % figRows = 2;
    % figCols = numel(ExtracellularVolumes);
    % figure; 
    % tiledlayout(figRows,figCols,'TileIndexing','columnmajor','TileSpacing','tight')
    % set(gcf,'color','w','position',[200 200 1600 950])
    % ylim_UB_WexAUC = zeros(numel(ExtracellularVolumes),1);
    % ylim_UB_WinAUC = zeros(numel(ExtracellularVolumes),1);
    % 
    % for i = 1:numel(ExtracellularVolumes)
    %     AUC1{i} = trapz(T{i},Y{i});
    %     ylim_UB_WexAUC(i,1) = max(AUC1{i}(Wex));
    %     ylim_UB_WinAUC(i,1) = max(AUC1{i}(Win));
    % end
    % 
    % for i = 1:numel(ExtracellularVolumes)
    %     nexttile 
    %     bar(AUC1{i}(Wex),'FaceColor','k','FaceAlpha',.1)
    %     hold on
    %     bar([ 1 ;nan],[AUC1{i}(Wex_deg) AUC1{i}(Wex_efflux); nan(1,2)],'stacked','FaceColor','flat');      
    %     hold off
    %     if i == 1
    %         ylabel('Extracellular Warhead (AUC)','fontsize',20,'FontWeight','bold')
    %     end
    %     if i > 1
    %         set(gca,'YTickLabel',[])
    %     end
    %     if i == numel(ExtracellularVolumes)
    %         h = legend({'',['Extracellular' newline 'Deconjugation'],'Efflux'},'Location','eastoutside');
    %         title(h,'Source')
    %     end
    %     ylim([0 max(ylim_UB_WexAUC)*11/10])
    %     set(gca,'FontSize',20,'XTickLabel',[])
    %     title(['10^-^' num2str(abs(log10(ExtracellularVolumes(i)))) ' L'],'fontsize',20,'FontWeight','bold')
    %     grid on; box on;
    % 
    %     nexttile
    %     bar(AUC1{i}(Win),'FaceColor','k','FaceAlpha',.1)
    %     hold on
    %     bar([ 1 ;nan],[AUC1{i}(Win_deg) AUC1{i}(Win_influx) AUC1{i}(Win_lys) AUC1{i}(Win_nuc); nan(1,4)],'stacked','FaceColor','flat');      
    %     hold off
    %     if i == 1
    %         ylabel('Intracellular Warhead (AUC)','fontsize',20,'FontWeight','bold')
    %     end
    %     if i > 1
    %         set(gca,'YTickLabel',[])
    %     end
    %     if i == numel(ExtracellularVolumes)
    %         h = legend({'',['Extracellular' newline 'Deconjugation'],'Influx','Lysosome','Nucleus'},'Location','eastoutside');
    %         title(h,'Source')
    %     end
    %     ylim([0 max(ylim_UB_WinAUC)*11/10])
    %     set(gca,'FontSize',20,'XTickLabel',[])
    %     grid on; box on;
    % end
    % 
    % drawnow;
    % fprintf("Done!\n")
end
    