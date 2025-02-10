doseType = "ADC"; % "ADC" or "Ab" or "PBD" or "IsotypeADC" or "None"
time = 24*7; % (hr) incubation time in hours
run(setup_file);

if plotFigure7B == "yes"

    fprintf("Generating Figure 7B...")

    % Simulate treated cells
    [T,Y] = ode23s(eqns_file,[0 time],InitCond,ode_options,p);

    % Stacked line plot
    figRows = 1;
    figCols = 2;

    fig7b = figure; % create figure
    tiledlayout(figRows,figCols,'TileSpacing','tight')
    nexttile
    area(T,Y(:,Wex),'linewidth',3,'FaceColor','k','FaceAlpha',.1)
    hold on
    area(T,Y(:,[Wex_deg, Wex_efflux]),'linewidth',3,'FaceAlpha',0.9)
    title('Extracellular Warhead','fontsize',20)
    ylabel('Concentration (nM)','fontsize',20,'FontWeight','bold')
    xlabel('Time (hours)','fontsize',20,'FontWeight','bold')
    xlim([0,time])
    h1 = legend({'',['Extracellular' newline 'Deconjugation'],'Efflux'},'Location','eastoutside');
    title(h1,'Source')
    set(gca,'FontSize',20)
    box on; grid on;

    nexttile
    area(T,Y(:,Win),'linewidth',3,'FaceColor','k','FaceAlpha',.1)
    hold on
    area(T,Y(:,[Win_deg, Win_influx, Win_lys, Win_nuc]),'linewidth',3,'FaceAlpha',0.9)
    title('Intracellular Warhead','fontsize',20)
    ylabel('Concentration (nM)','fontsize',20,'FontWeight','bold')
    xlabel('Time (hours)','fontsize',20,'FontWeight','bold')
    xlim([0,time])
    h2 = legend({'',['Extracellular' newline 'Deconjugation'],'Influx','Lysosome','Nucleus'},'Location','eastoutside');
    title(h2,'Source')
    set(gca,'FontSize',20)
    box on; grid on;
    set(gcf,'color','w','position',[200 120 1300 400])
    
    drawnow;
    fprintf("Done!\n")

end

if plotFigure7C == "yes"

    fprintf("Generating Figure 7C...")

    p.kdeconjugation = kdeconjugation_MouseSerum;

    % Simulate treated cells
    [T,Y] = ode23s(eqns_file,[0 time],InitCond,ode_options,p);

    % Stacked line plot
    figRows = 1;
    figCols = 2;

    fig7c = figure; % create figure
    tiledlayout(figRows,figCols,'TileSpacing','tight')

    nexttile
    area(T,Y(:,Wex),'linewidth',3,'FaceColor','k','FaceAlpha',.1)
    hold on
    area(T,Y(:,[Wex_deg, Wex_efflux]),'linewidth',3,'FaceAlpha',0.9)
    title('Extracellular Warhead','fontsize',20)
    ylabel('Concentration (nM)','fontsize',20,'FontWeight','bold')
    xlabel('Time (hours)','fontsize',20,'FontWeight','bold')
    xlim([0,time])
    h1 = legend({'',['Extracellular' newline 'Deconjugation'],'Efflux'},'Location','eastoutside');
    title(h1,'Source')
    set(gca,'FontSize',20)
    box on
    grid on

    nexttile
    area(T,Y(:,Win),'linewidth',3,'FaceColor','k','FaceAlpha',.1)
    hold on
    area(T,Y(:,[Win_deg, Win_influx, Win_lys, Win_nuc]),'linewidth',3,'FaceAlpha',0.9)
    title('Intracellular Warhead','fontsize',20)
    ylabel('Concentration (nM)','fontsize',20,'FontWeight','bold')
    xlabel('Time (hours)','fontsize',20,'FontWeight','bold')
    xlim([0,time])
    h2 = legend({'',['Extracellular' newline 'Deconjugation'],'Influx','Lysosome','Nucleus'},'Location','eastoutside');
    title(h2,'Source')
    set(gca,'FontSize',20)
    box on
    grid on
    set(gcf,'color','w','position',[200 120 1300 400])

    drawnow;
    fprintf("Done!\n")

end