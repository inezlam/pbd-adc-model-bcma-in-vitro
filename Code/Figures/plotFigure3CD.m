doseType = "ADC"; % "ADC" or "Ab" or "PBD" or "IsotypeADC" or "None"
run(setup_file);

if plotFigure3C == "yes"

    fprintf("Generating Figure 3C...")

    % Experimental conditions
    time = 24*4; % (hr) incubation time in hours
    Vwell = 80e-6; % volume of media added per well for cells
    Vdose = 20e-6; % volume of PBD dose
    V = Vwell + Vdose; % Volume per well in cytotoxicity assay plus volume of PBD added per well % (L) Total Volume of compartment
    Ncell = 5000; % # cells in 2 mL solution
    p.Vmedia = V - (Vcell*Ncell); % (L) working volume of media in 96 well plate (taken for expt)
    InitCond(Cells) = Ncell;
    ADCs = {"MEDI2228"};
    ADC_doses = [1.85999E-05 5.58656E-05 0.000167531 0.000502644 0.001507935 0.004523797 0.013571392 0.040714182 0.122142546 0.366427639 1.099283048 3.297848484]; % nM - initial concentrations of ADC        
    x_exp = ADC_doses;
    y_exp = [92.33 103.55 97.36 106.42 98.37 102.4 70.33 50.05 29.77 16.68 6.62 2.44];
    SD = [1.22 6.92 10.37 10.17 6.1 26.44 4.27 6.51 4.27 3.25 2.03 0.2];
    
    % Run Simulation
    numSimPts = 12*5;
    x_sim = logspace(log10(ADC_doses(1,1)),log10(ADC_doses(1,end)),numSimPts); % nM
    [~,y_sim] = dose_response(eqns_file, p, time, x_sim, ADC, InitCond);
    
     % Plot figure
    figure;
    hold on;
    plot(log10(x_sim),y_sim,'linewidth',4,'color','k','DisplayName','Simulated')
    errorbar(log10(x_exp),y_exp,SD,'o',"MarkerSize",10,'Color',"#5DA9D8",'MarkerEdgeColor',"#5DA9D8","MarkerFaceColor",'#5DA9D8')
    ylabel('Cell Survival (%)','fontsize',20,'FontWeight','bold')
    xlabel('log_{10}[ADC] (nM)','fontsize',20,'FontWeight','bold')
    set(gca,'FontSize',20)
    set(gcf,'color','w','position',[200 120 500 400])
    box on
    grid on
    xlim([-5 1])
    ylim([0 130])
    [~,objh] = legend({'Simulated','Experimental'},'Location','southwest','FontSize',18);
    drawnow;
    fprintf("Done!\n")

end

if plotFigure3D == "yes"

    fprintf("Generating Figure 3D...")

    % Experimental conditions
    time = 24*4; 
    sBCMA_ng_mL = [0 75 270 720];
    Ag_s0_multiplier = (1e-6 / Ag_s_MW) * 10^9;
    sBCMA_nM = sBCMA_ng_mL.*Ag_s0_multiplier;

    V = 100e-6; % Volume per well in cytotoxicity assay plus volume of PBD added per well % (L) Total Volume of compartment
    Ncell = 5000; % # cells in 2 mL solution
    p.Vmedia = V - (Vcell*Ncell); % (L) working volume of media in 96 well plate (taken for expt)
    InitCond(Cells) = Ncell;

    ADC_doses = [1.85999E-05 5.58656E-05 0.000167531 0.000502592 0.001510415 0.004524648 0.013573944 0.040715237 0.122145712 0.366430541 1.099285026 3.297848484];
    x_exp = ADC_doses;
    y_exp = [101.97 101.24 96.09 106.32 98.11 93.7 56.28 38.52 23.27 14.58 7.84 3.61;  % 0 ng/mL
            99.6 110.79 102.82 97.64 88.16 87.72 60.92 40.83 26.04 17 9.23 4.17; % 75 ng/mL
            95.55 113.79 95.39 90.9 92.04 89.28 64.39 41.77 23.75 16.02 9.47 4.71; % 270 ng/mL
            106 114.61 102.43 103.13 101.74 97.48 67.13 46.98 27.28 16.88 10.15 5.5]; % 720 ng/mL
    SD = [9.09 1.99 4.59 16.11 7.45 6.41 4.07 0.09 2.08 1.73 0.17 0.09; % 0 ng/mL
        2.86 1.16 3.31 8.67 4.38 5.72 1.97 4.83 0.89 1.16 0.36 0.89; % 75 ng/mL
        9.34 10.48 4.67 6.12 4.67 2.91 1.07 1.38 0.99 1.53 0.08 0.38; % 270 ng/mL
        2.31 3.57 3.01 6.09 1.75 5.81 0.56 2.73 1.33 0.63 0.21 1.19]; % 720 ng/mL

    numSimPts = 8*8;

    x = cell(length(sBCMA_nM),1);
    y = cell(length(sBCMA_nM),1);
    x_sim = logspace(log10(ADC_doses(1,1)),log10(ADC_doses(1,end)),numSimPts); 

    for i = 1:length(sBCMA_nM)
        InitCond(Ag_s) = sBCMA_nM(i);
        [x{i},y{i}] = dose_response(eqns_file,p, time, x_sim, ADC, InitCond);
    end

    colors = ["#A1E5E1","#86BFBC","#597F7D","#364C4C"];
    figure
    hold on
    for i = 1:length(sBCMA_nM)
        color = colors(i);
        plot(log10(x{i}),y{i},'linewidth',2.5,'DisplayName',['Sim: ',num2str(sBCMA_ng_mL(i)),' ng/mL'],'Color',color)
        errorbar(log10(x_exp),y_exp(i,:),SD(i,:),'o',"MarkerSize",8,...
            'DisplayName',['Exp: ',num2str(sBCMA_ng_mL(i)),' ng/mL'],...
            'Color',color,'MarkerEdgeColor',color,"MarkerFaceColor",color)
    end
    set(gca,'FontSize',20)
    ylabel('Cell Survival (%)','fontsize',20,'FontWeight','bold')
    xlabel('log_{10}[ADC] (nM)','fontsize',20,'FontWeight','bold')
    set(gcf,'color','w','position',[200 120 500 400])
    box on
    grid on
    xlim([-5 1])
    ylim([0 130])

    % Set legend
    colorNames = {'0 ng/mL','75 ng/mL','270 ng/mL','720 ng/mL'};
    hold on
    h = [];
    for j = numel(colors):-1:1
        h(j) = plot(NaN',NaN','Color',colors{j},'LineStyle','-','LineWidth',2.5,'Marker','o','MarkerFaceColor',colors(j),'MarkerEdgeColor',colors(j),'MarkerSize',8);
    end
    
    hold off
    leg = legend(h,colorNames,'Location','southwest','FontSize',18);
    title(leg,'sBCMA Levels')
    drawnow;
    fprintf("Done!\n")
    
end