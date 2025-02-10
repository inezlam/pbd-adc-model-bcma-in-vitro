doseType = "ADC";
time = 24*4; % (hr) incubation time in hours
run(setup_file);

if plotFigure5A == "yes"

    fprintf("Generating Figure 5A...\n")

    %% DAR Simulations
    DAR_values = [1 2 3 4 5 6 7 8];
    
    figure;
    hold on;
    cm = flip(colormap(winter(size(DAR_values',1))));
    
    for i = 1:numel(DAR_values)
        
        p_test = p;
        p_test.DAR = DAR_values(i);
    
        % Simulate treated cells
        % fprintf("Simulate treated cells... ")
        [T1,Y1] = ode23s(eqns_file,[0 time],InitCond,ode_options,p_test);
        % fprintf("Done!\n")
        
        % Simulate untreated cells (control)
        % fprintf("Simulate untreated cells... ")
        [T_control,Y_control] = ode23s(eqns_file,T1,InitCond_noADC,ode_options,p_test);        
        CellPop_control = Y_control(:,length(InitCond));   
        % fprintf("Done!\n") 
    
        plot(T1,(Y1(:,Cells)./CellPop_control)*100,'linewidth',2,'Color', cm(i,:),'DisplayName',num2str(DAR_values(i)))
        disp(['Finished loop ',num2str(i),' out of ',num2str(numel(DAR_values))])   
        
    end
            
    ylabel('Cell Survival (%)','fontsize',20,'FontWeight','bold')
    xlabel('Time (hours)','fontsize',20,'FontWeight','bold')
    xlim([0,time])
    ylim([0 100])
    set(gca,'FontSize',20)
    set(gcf,'color','w','position',[200 120 500 400])
    box on
    grid on
    leg = legend('show','Location','southwest','FontSize',18);
    title(leg,'DAR')
    
    fprintf("Done!\n")
    drawnow;

end

if plotFigure5B == "yes"

    fprintf("Generating Figure 5B...\n")

    %% DAR Simulations
    DAR_values = [1 2 3 4 5 6 7 8];
    
    figure;
    hold on;
    cm = flip(colormap(winter(size(DAR_values',1))));
    
    for i = 1:numel(DAR_values)
        
        p_test = p;
        p_test.DAR = DAR_values(i);
    
        % Simulate treated cells
        [T1,Y1] = ode23s(eqns_file,[0 time],InitCond,ode_options,p_test);
        plot(T1,(Y1(:,Win_influx)),'linewidth',2,'Color', cm(i,:),'DisplayName',num2str(DAR_values(i))) 
        disp(['Finished loop ',num2str(i),' out of ',num2str(numel(DAR_values))])   
        
    end
    
    % title('Parameter Scan for DAR','fontsize',20)
    
    ylabel('Influxed Warhead (nM)','fontsize',20,'FontWeight','bold')
    xlabel('Time (hours)','fontsize',20,'FontWeight','bold')
    xlim([0,time])
    % ylim([0 100])
    set(gca,'FontSize',20)
    set(gcf,'color','w','position',[200 120 500 400])
    box on
    grid on
    leg = legend('show','Location','northwest','FontSize',18);
    title(leg,'DAR')
    
    fprintf("Done!\n")
    drawnow;
end

if plotFigure5C == "yes"

    fprintf("Generating Figure 5C...\n")

    %% DAR Simulations
    DAR_values = [1 2 3 4 5 6 7 8];
    numSimPts = 8*8;
    ADC_doses = logspace(-5,1,numSimPts); %*ADC0/p.Vmedia;
    y0 = InitCond; % nM
    x = cell(numel(DAR_values),1);
    y = cell(numel(DAR_values),1);
    
    figure;
    hold on;
    cm = flip(colormap(winter(size(DAR_values',1))));
    
    for i = 1:numel(DAR_values)
        
        p_test = p;
        p_test.DAR = DAR_values(i);
        [x{i},y{i}] = dose_response(eqns_file, p_test, time, ADC_doses, ADC, y0);
        plot(log10(x{i}),y{i},'Color', cm(i,:),'linewidth',2,'DisplayName',[num2str(DAR_values(i))])        
        disp(['Finished loop ',num2str(i),' out of ',num2str(numel(DAR_values))])   
        
    end

    ylabel('Cell Survival (%)','fontsize',20,'FontWeight','bold')
    xlabel('log_{10}[ADC] (nM)','fontsize',20,'FontWeight','bold')
    set(gca,'FontSize',20)
    set(gcf,'color','w','position',[200 120 500 400])
    box on
    grid on
    xlim([-5 0.5])
    leg = legend('show');
    title(leg,'DAR')

    fprintf("Done!\n")
    drawnow;

end

if plotFigure5D == "yes"

    fprintf("Generating Figure 5D...\n")

    ADC_doses = logspace(-5,1,7);

    % Define DAR values
    DAR_values = linspace(1,8,8);
    Y = zeros(numel(DAR_values),length(ADC_doses));

    % For each DAR value
    for i = 1:numel(DAR_values)

        p_test = p;
        p_test.DAR = DAR_values(i);

        % Run simulation - cell survival at X time, X dose
        [~,Y(i,:)] = dose_response(eqns_file, p_test, time, ADC_doses, ADC, InitCond);

        disp(['Finished loop ',num2str(i),' out of ',num2str(numel(DAR_values))])   

    end

    % Define colormap
    cm_r = linspace(0.8,0/255,length(ADC_doses));
    cm_g = linspace(0.9,50/255,length(ADC_doses));
    cm_b = linspace(1,0.8,length(ADC_doses));
    cm = [cm_r' cm_g' cm_b'];

    % Plot figure
    figure;
    hold on;

    for i = 1:numel(ADC_doses)
        plot(DAR_values,Y(:,i),'linewidth',2.5,'LineStyle','-','Marker','o','Color',cm(i,:),'MarkerFaceColor',cm(i,:),'MarkerEdgeColor',cm(i,:),'MarkerSize',7)
    end
    ylabel('Cell Survival (%)','fontsize',20,'FontWeight','bold')
    xlabel('DAR','fontsize',20,'FontWeight','bold')
    set(gca,'FontSize',20)
    set(gcf,'color','w','position',[200 120 600 400])
    box on
    grid on
    xlim([min(DAR_values) max(DAR_values)])
    leg = legend({'10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}','10^{1}'},'Location','eastoutside','FontSize',18);
    title(leg,'Dose (nM)')
    
    fprintf("Done!\n")
    drawnow;

end

