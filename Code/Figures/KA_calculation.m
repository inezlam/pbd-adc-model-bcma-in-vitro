function [KA,T0,Y0] = KA_calculation(fxn,p,ExposureTime,IncubationTime,dose,InitCond,DrugSpecies,CrosslinkSpecies,plotCrosslinkFigure)

%     fprintf("Calculating value of KA (half maximal concentration of Wn-DNA for crosslinking)... ")
    InitCond(DrugSpecies) = dose; % nM, calculated value for Wex to achieve 50% crosslinking
    [T0,Y0] = exposure_incubation(fxn, p, ExposureTime, IncubationTime, InitCond, DrugSpecies);
    KA = Y0(end,CrosslinkSpecies);
    
    %% Check function
    if plotCrosslinkFigure == "Yes"
        Crosslinks = 1 ./ (1 + (KA./Y0(:,CrosslinkSpecies)).^p.n);
        figure
        plot(T0,Crosslinks,'linewidth',2)
        hold on
        xline(2,'k-')
        yline(0.5,'r--','linewidth',2)
        title(['Crosslink Formation: KA = ' num2str(KA)],'fontsize',20)
        ylabel('Fraction Crosslinked','fontsize',20)
        xlabel('Time (hours)','fontsize',20)
        xlim([0,IncubationTime+ExposureTime])
        set(gca,'FontSize',20)
        set(gcf,'color','w')
        grid on; box on;
    end
    
%     fprintf("The calculated value for KA is %e nM.\n",KA)
%     toc
    
end  