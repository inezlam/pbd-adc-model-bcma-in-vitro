doseType = "IsotypeADC"; % "ADC" or "Ab" or "PBD" or "IsotypeADC" or "None"
run(setup_file);

if plotFigure3A == "yes"

    fprintf("Generating Figure 3A...")

    % Experimental Conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    time = 96; % (hr) incubation time in hours
    Vwell = 80e-6; % volume of media added per well for cells
    Vdose = 20e-6; % volume of PBD dose
    V = Vwell + Vdose; % Volume per well in cytotoxicity assay plus volume of PBD added per well % (L) Total Volume of compartment
    Ncell = 5000; % # cells in 2 mL solution
    p.Vmedia = V - (Vcell*Ncell); % (L) working volume of media in 96 well plate (taken for expt)
    InitCond(Cells) = Ncell;
    % Experimental Data Points - combine and calculate average for two experiments
    ADC_doses = [1.85999E-05 5.58656E-05 0.000167531 0.000502644 0.001507935 0.004523797 0.013571392 0.040714182 0.122142546 0.366427639 1.099283048 3.297848484]; % nM - initial concentrations of ADC
    x_exp = ADC_doses;
    y_exp = [106.195 112.21 101.33 109.155 96.435 105.645 105.44 107.305 98.83 108.025 93.49 75.445]; % Avg of 2 expts at 96 hours
    y_exp1 = [106.85 105.85 94.37 104.52 86.55 105.69 110.85 108.35 101.36 116.67 99.53 85.71];
    SD1 = [1.88 2.35 0.71 4.71 3.77 11.06 2.82 20.95 2.12 23.77 4.24 2.12];
    n1 = 2;
    y_exp2 = [105.54 118.57 108.29 113.79 106.32 105.6 100.03 106.26 96.3 99.38 87.45 65.18];
    SD2 = [8.06 5.74 7.69 8.43 14.36 6.86 5.84 3.34 0.74 0.28 3.43 0.28];
    n2 = 2;
    var1 = SD1.^2;
    var2 = SD2.^2;
    q1 = (n1-1).*var1 + n1.*y_exp1.^2;
    q2 = (n2-1).*var2 + n2.*y_exp2.^2;
    qc = q1 + q2;
    SD = sqrt((qc - (n1+n2).*y_exp.^2)./(n1+n2-1));
    SDc = sqrt(var1 + var2);

    % Run simulation
    numSimPts = 12*5;
    x_sim = logspace(log10(ADC_doses(1,1)),log10(ADC_doses(1,end)),numSimPts); % nM
    [~,y_sim] = dose_response(eqns_file, p, time, x_sim, ADC, InitCond);
    
    % Plot figure
    figure;
    set(gcf,'color','w','position',[200 120 500 400])
    hold on;
    plot(log10(x_sim),y_sim,'linewidth',4,'color','k','DisplayName','Simulated')
    errorbar(log10(x_exp),y_exp,SD,'o',"MarkerSize",10,'Color',"#5DD9C2",'MarkerEdgeColor',"#5DD9C2","MarkerFaceColor",'#5DD9C2')
    ylabel('Cell Survival (%)','fontsize',20,'FontWeight','bold')
    xlabel('log_{10}[Isotype Control ADC] (nM)','fontsize',20,'FontWeight','bold')
    set(gca,'FontSize',20)
    box on
    grid on
    xlim([-5 1])
    ylim([0 130])
    [~,objh] = legend({'Simulated','Experimental'},'Location','southwest','FontSize',18);
    drawnow;
    fprintf("Done!\n")
end

if plotFigure3B == "yes"

    fprintf("Generating Figure 3B...")

    % Experimental Conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    time = 72; % (hr) incubation time in hours
    Vwell = 80e-6; % volume of media added per well for cells
    Vdose = 20e-6; % volume of PBD dose
    V = Vwell + Vdose; % Volume per well in cytotoxicity assay plus volume of PBD added per well % (L) Total Volume of compartment
    Ncell = 5000; % # cells in 2 mL solution
    p.Vmedia = V - (Vcell*Ncell); % (L) working volume of media in 96 well plate (taken for expt)
    InitCond(Cells) = Ncell;
    ADC_doses = [0.000501273 0.001501701 0.00455103 0.0135871 0.04069545 0.122142546 0.366427639 1.099283048 3.297848484]; % nM - initial concentrations of ADC
    x_exp = ADC_doses;
    y_exp = [108.1 96.2 92.3 105.3 115 98.7 96.3 89.2 91.2];
    SD = [24.2 7 17.1 12.1 2.4 3.9 11.9 0.7 3.6];
    numSimPts = 20;
    x_sim = logspace(log10(ADC_doses(1,1)),log10(ADC_doses(1,end)),numSimPts); % nM
    [~,y_sim] = dose_response(eqns_file, p, time, x_sim, ADC, InitCond);

    % Plot figure
    figure;
    hold on;
    plot(log10(x_sim),y_sim,'linewidth',4,'color','k','DisplayName','Simulated')
    errorbar(log10(x_exp),y_exp,SD,'o',"MarkerSize",10,'Color',"#5DD9C2",'MarkerEdgeColor',"#5DD9C2","MarkerFaceColor",'#5DD9C2')
    ylabel('Cell Survival (%)','fontsize',20,'FontWeight','bold')
    xlabel('log_{10}[Isotype Control ADC] (nM)','fontsize',20,'FontWeight','bold')
    set(gca,'FontSize',20)
    set(gcf,'color','w','position',[200 120 500 400])
    box on
    grid on
    xlim([-3.5 1])
    ylim([0 135])
    [~,objh] = legend({'Simulated','Experimental'},'Location','southwest','FontSize',18);
    drawnow;
    fprintf("Done!\n")

end