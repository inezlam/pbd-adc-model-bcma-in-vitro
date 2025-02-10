doseType = "PBD";
run(setup_file);

if plotFigure2B == "yes" % Optimization for KA for PBD - crosslinks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf("Generating Figure 2B...")

    % Inital Conditions - No ADC, starting with extracellular warhead only
    InitCond(ADC) = 0;
    x_exp = [0.17139573 0.51286138 1.68267406 5.08159443 16.7494288];
    PBDs = {"SG3199"};
    x = cell(numel(PBDs),1);
    y = cell(numel(PBDs),1);
    optimal_KA = nan(numel(PBDs),1);       
    y_exp = [22.383748 31.6885669 49.6900787 81.7669291 94.8552441]./100;
    SD = [12.52677921 17.31168629 2.473098042 3.924944559 7.849398133]./100;
    ExposureTime = 2;
    IncubationTime = 24;
    time = ExposureTime + IncubationTime;

    %% Print simulation figure with more points
    p_test = p;
    numSimPts = 5*8;
    x_sim = logspace(-1,1.5,numSimPts); % nM
    [~,y_sim] = dose_response_crosslink_incubation(eqns_file, p_test, ExposureTime, IncubationTime, x_sim, Wex, Wn_DNA, InitCond_PBD); % with incubation
    
    % Plot final figure
    color = "#EDB210";
    figure;
    hold on;
    plot(log10(x_sim),y_sim,'linewidth',4,'color','k','DisplayName','Simulated');
    errorbar(log10(x_exp),y_exp,SD,'o',"MarkerSize",10,'DisplayName','Experimental',...
        'Color',color,'MarkerEdgeColor',color,"MarkerFaceColor",color)
    ylabel('Fraction Crosslinked','fontsize',20,'FontWeight','bold')
    xlabel('log_{10}[PBD] (nM)','fontsize',20,'FontWeight','bold')
    set(gca,'FontSize',20)
    set(gcf,'color','w','position',[100 100 500 400])
    xlim([-1 1.5])
    box on
    grid on
    [~,objh] = legend({'Simulated','Experimental'},'Location','northwest','FontSize',18); % Set legend marker size for scatter points
    drawnow;
    fprintf("Done!\n")
    
end

if plotFigure2C == "yes" % Crosslink Formation over Time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf("Generating Figure 2C...")

    % Experimental conditions
    ExposureTime = 2;
    IncubationTime = 36;
    time = ExposureTime + IncubationTime;
    InitCond = InitCond_noADC;
    InitCond(Wex) = 1.69; % Dose given in Hartley 2018 paper
    t_exp = [0 2 4 8 12 24 36] + 2;
    y_exp = [43.5227273 40.7954545 53.75 56.1363636 65 65 61.5909091]./100;
    SD = [18.13807531 29.18410042 13.11715481 17.0083682 11.79916318 13.36820084 17.57322176]./100;
  
    % Simulate exposure and iincubation
    [T,Y] = exposure_incubation(eqns_file, p, ExposureTime, IncubationTime, InitCond, Wex);
    Crosslinks = 1 ./ (1 + (p.KA./Y(:,Wn_DNA)).^p.n);
           
    % Plot final figure
    color = "#EDB210";
    figure;
    hold on;
    xl = xline(2,'--r','Washout','FontSize',18,'LineWidth',2);
    plot(T,Crosslinks,'linewidth',4,'color','k','DisplayName','Simulated');
    eb = errorbar(t_exp,y_exp,SD,'o',"MarkerSize",10,'DisplayName','Experimental',...
        'Color',color,'MarkerEdgeColor',color,"MarkerFaceColor",color);
    ylabel('Fraction Crosslinked','fontsize',20,'FontWeight','bold')
    xlabel('Time (hours)','fontsize',20,'FontWeight','bold')
    set(gca,'FontSize',20)
    set(gcf,'color','w','position',[200 120 500 400])
    xlim([0 38])
    ylim([0 1])
    box on
    grid on
    [~,objh] = legend({'','Simulated','Experimental'},'Location','northeast','FontSize',18); % Set legend marker size for scatter points
    uistack(eb,'top')
    drawnow;
    fprintf("Done!\n")

end

if plotFigure2D == "yes"

    fprintf("Generating Figure 2D...")

    % Experimental conditions
    Vwell = 80e-6; % volume of media added per well for cells
    Vdose = 20e-6; % volume of PBD dose
    V = Vwell + Vdose; % Volume per well in cytotoxicity assay plus volume of PBD added per well       
    cellDensity = 5000/80e-6; % cells/mL (from PBD cytotoxicity assay)
    Ncell = 5000; % # cells in 2 mL solution
    p.Ncell = Ncell;
    p.Vmedia = V - (Vcell*Ncell); % (L) working volume of media in 96 well plate (taken for expt)
    Vmedia = p.Vmedia;
    PBD_doses = [7.8125E-12 1.5625E-10 3.1250E-09 6.2500E-08 1.2500E-06 2.5000E-05 5.0000E-04 1.0000E-02]./V; % nM - initial concentrations of Wex 
    x_exp = PBD_doses;
    y_exp = [102.938 107.498 108.493 104.293 30.094 9.81 2.819 3.786];
    SD = [1.993 1.172 2.267 0.703 0.586 0.352 0.938 0.664];

    % Run Simulations
    time = 24*3; % (hr) incubation time in hours
    numSimPts = 8*5;
    InitCond(ADC) = 0; % Inital Conditions - No ADC, starting with extracellular warhead only
    PBD_doses = logspace(-8,3,numSimPts);
    PBD_doses_plot = logspace(-7,2,numSimPts); %*PBD0/Vwell; % nM
    [x,y] = dose_response(eqns_file, p, time, PBD_doses, Wex, InitCond);

    % Plot figure
    color = "#7030A0";
    figure;
    hold on;
    plot(log10(PBD_doses),y,'linewidth',4,'color','k','DisplayName','Simulated')
    errorbar(log10(x_exp),y_exp,SD,'o',"MarkerSize",10,'DisplayName','Experimental',...
        'Color',color,'MarkerEdgeColor',color,"MarkerFaceColor",color)
    ylabel('Cell Survival (%)','fontsize',20,'FontWeight','bold')
    xlabel('log_{10}[PBD] (nM)','fontsize',20,'FontWeight','bold')
    set(gca,'FontSize',20)
    set(gcf,'color','w','position',[200 120 500 400])
    box on
    grid on
    xlim([-8 3])
    [~,objh] = legend({'Simulated','Experimental'},'Location','northeast','FontSize',18); % Set legend marker size for scatter points
    drawnow;
    fprintf("Done!\n")
end
