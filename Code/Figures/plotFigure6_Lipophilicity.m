doseType = "ADC"; % "ADC" or "Ab" or "PBD" or "IsotypeADC" or "None"
time = 24*4; % (hr) incubation time in hours
run(setup_file);

%% Varying R, Cell Survival vs. Time
if plotFigure6C == "yes" 

    fprintf("Generating Figure 6C...\n")

    R_values = logspace(-1,3,5);
    figure;
    hold on;
    cm = flip(colormap(winter(size(R_values',1))));
    
    for i = 1:numel(R_values)
        
        p_test = p;
        p_test.keff = p_test.kinf/(1 + R_values(i)); % (hr-1) efflux rate constant
    
        % Simulate treated cells
        [T1,Y1] = ode23s(eqns_file,[0 time],InitCond,ode_options,p_test);
        
        % Simulate untreated cells (control)
        [T_control,Y_control] = ode23s(eqns_file,T1,InitCond_noADC,ode_options,p_test);        
        CellPop_control = Y_control(:,length(InitCond));   
    
        plot(T1,(Y1(:,Cells)./CellPop_control)*100,'linewidth',2,'Color', cm(i,:),'DisplayName',num2str(R_values(i)))        
        disp(['Finished loop ',num2str(i),' out of ',num2str(numel(R_values))])   
        
    end
    % Label figure
    ylabel('Cell Survival (%)','fontsize',20,'FontWeight','bold')
    xlabel('Time (hours)','fontsize',20,'FontWeight','bold')
    xlim([0,time])
    ylim([0 100])
    set(gca,'FontSize',20)
    set(gcf,'color','w','position',[200 120 500 400])
    box on
    grid on
    leg = legend('show','Location','southwest','FontSize',18);
    title(leg,'R')
    drawnow;
    fprintf("Done!\n")

end

if plotFigure6D == "yes"

    fprintf("Generating Figure 6D...\n")

    R_values = logspace(-1,3,5);
    numSimPts = 8*8;
    ADC_doses = logspace(-5,1,numSimPts);
    y0 = InitCond; % nM
    x = cell(numel(R_values),1);
    y = cell(numel(R_values),1);        
    
    figure;
    hold on;
    cm = flip(colormap(winter(size(R_values',1))));
    
    for i = 1:numel(R_values)
        
        p_test = p;
        p_test.keff = p_test.kinf/(1 + R_values(i)); % (hr-1) efflux rate constant
        [x{i},y{i}] = dose_response(eqns_file, p_test, time, ADC_doses, ADC, y0);
        plot(log10(x{i}),y{i},'Color', cm(i,:),'linewidth',2,'DisplayName',[num2str(R_values(i))])      
        disp(['Finished loop ',num2str(i),' out of ',num2str(numel(R_values))])   
        
    end
    
    ylabel('Cell Survival (%)','fontsize',20,'FontWeight','bold')
    xlabel('log_{10}[ADC] (nM)','fontsize',20,'FontWeight','bold')
    set(gca,'FontSize',20)
    set(gcf,'color','w','position',[200 120 500 400])
    box on
    grid on
    xlim([-5 0.5])
    leg = legend('show','Location','southwest','FontSize',18);
    title(leg,'R')
    drawnow;
    fprintf("Done!\n")
    
end

%% Varying R, Influxed Warhead vs. Time
if plotFigure6E == "yes"

    fprintf("Generating Figure 6E...\n")
    
    R_values = logspace(-1,3,5);
    
    figure;
    hold on;
    cm = flip(colormap(winter(size(R_values',1))));
    
    for i = 1:numel(R_values)
        
        p_test = p;
        p_test.keff = p_test.kinf/(1 + R_values(i)); % (hr-1) efflux rate constant
    
        % Simulate treated cells
        [T1,Y1] = ode23s(eqns_file,[0 time],InitCond,ode_options,p_test);
        plot(T1,(Y1(:,Win_influx)),'linewidth',2,'Color', cm(i,:),'DisplayName',num2str(R_values(i)))        
        disp(['Finished loop ',num2str(i),' out of ',num2str(numel(R_values))])   
        
    end
    
    ylabel('Influxed Warhead (nM)','fontsize',20,'FontWeight','bold')
    xlabel('Time (hours)','fontsize',20,'FontWeight','bold')
    xlim([0,time])
    set(gca,'FontSize',20)
    set(gcf,'color','w','position',[200 120 500 400])
    box on
    grid on
    leg = legend('show','Location','northwest','FontSize',18);
    title(leg,'R')
    drawnow;
    fprintf("Done!\n")

end