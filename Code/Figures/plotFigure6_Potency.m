doseType = "ADC"; 
time = 24*4; % (hr) incubation time in hours
run(setup_file);

if plotFigure6A == "yes"

    fprintf("Generating Figure 6A...\n")

    %% kkill Simulations
    kkill_values = [1e-3 1e-2 1e-1 1e0 1e1 1e2];
    
    figure;
    hold on;
    cm = flip(colormap(winter(size(kkill_values',1))));
    
    for i = 1:numel(kkill_values)
        
        p_test = p;
        p_test.kkill = kkill_values(i);
        
        % Simulate treated cells
        [T1,Y1] = ode23s(eqns_file,[0 time],InitCond,ode_options,p_test);
        
        % Simulate untreated cells (control)
        [T_control,Y_control] = ode23s(eqns_file,T1,InitCond_noADC,ode_options,p_test);        
        CellPop_control = Y_control(:,length(InitCond));   
    
        plot(T1,(Y1(:,Cells)./CellPop_control)*100,'linewidth',2,'Color', cm(i,:),'DisplayName',num2str(kkill_values(i)))
        
        disp(['Finished loop ',num2str(i),' out of ',num2str(numel(kkill_values))])   
        
    end
    
    ylabel('Cell Survival (%)','fontsize',20,'FontWeight','bold')
    xlabel('Time (hours)','fontsize',20,'FontWeight','bold')
    xlim([0,time])
    ylim([0 100])
    set(gca,'FontSize',20)
    set(gcf,'color','w','position',[200 120 500 400])
    box on
    grid on
    leg = legend('show','Location','southeast','FontSize',18);
    title(leg,'k_k_i_l_l')
    drawnow;
    fprintf("Done!\n")

end

if plotFigure6B == "yes"

    fprintf("Generating Figure 6B...\n")

    %% Simulations
    kkill_values = [1e-3 1e-2 1e-1 1e0 1e1 1e2];
    numSimPts = 8*8;
    ADC_doses = logspace(-5,1,numSimPts); %*ADC0/p.Vmedia;
    y0 = InitCond; % nM
    x = cell(numel(kkill_values),1);
    y = cell(numel(kkill_values),1);

    figure;
    hold on;
    cm = flip(colormap(winter(size(kkill_values',1))));
    
    for i = 1:numel(kkill_values)
        
        p_test = p;
        p_test.kkill = kkill_values(i);
        [x{i},y{i}] = dose_response(eqns_file, p_test, time, ADC_doses, ADC, y0);
        plot(log10(x{i}),y{i},'Color', cm(i,:),'linewidth',2,'DisplayName',num2str(kkill_values(i)))        
        disp(['Finished loop ',num2str(i),' out of ',num2str(numel(kkill_values))])   
        
    end
    
    ylabel('Cell Survival (%)','fontsize',20,'FontWeight','bold')
    xlabel('log_{10}[ADC] (nM)','fontsize',20,'FontWeight','bold')
    set(gca,'FontSize',20)
    set(gcf,'color','w','position',[500 500 500 400])
    xlim([-5 0.5])
    leg = legend('show','Location','southeast','FontSize',18);
    title(leg,'k_k_i_l_l')
    grid on
    box on
    drawnow;
    fprintf("Done!\n")

end