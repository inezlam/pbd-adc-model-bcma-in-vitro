%% Optimization for Ag Production Rate
tic
fprintf('Starting optimization for Ag production rate...\n')

ode_options = odeset('MaxStep',5e-1, 'AbsTol', 1e-12,'RelTol', 1e-5,'InitialStep', 1e-2);

time_ss = 24*7*2*2; % (hr) incubation time in hours
numTimePoints = 50;
Ag_ss = InitCond(Ag);
[~, ~, Ag_endo_ss, Ag_lys_ss] = Ag_ss_calculation2(eqns_file,ode_options,p,time_ss,Ag,Ag_endo,Ag_lys,InitCond_noADC,p.kAgProduction);
InitCond_noADC(Ag_endo) = Ag_endo_ss(end);
InitCond_noADC(Ag_lys) = Ag_lys_ss(end);

% Optimize 1 param: kAgProduction
p_test = p;
p_test.kg = 0; % Set cell growth rate to 0 to find Ag production rate at static cell population
AgProd_0 = p_test.kAgProduction; % initial guess
lb = 1e-20;
ub = 1e-7; 
opt_options = optimoptions('lsqnonlin','Display','iter');   
opt_options.OptimalityTolerance = 1.000000e-50;
opt_options.FunctionTolerance = 1.000000e-50;
opt_options.StepTolerance = 1.000000e-50;
cost = @(AgProd_0) costfxn_AgProd(eqns_file, ode_options, AgProd_0, Ag_ss, p_test, time_ss, Ag, Ag_endo, Ag_lys, InitCond_noADC, numTimePoints);
optimal_AgProduction = lsqnonlin(cost,AgProd_0,lb,ub,opt_options);
% fprintf("The 1st optimized value for kAgProduction is %e nmol/hr.\n",optimal_AgProduction)

% Calculate Ag_endo_ss
[T_ss, Ag_ss_new, Ag_endo_ss_new, Ag_lys_ss_new] = Ag_ss_calculation(eqns_file,ode_options,p_test,time_ss,Ag,Ag_endo,Ag_lys,InitCond_noADC,optimal_AgProduction);
% InitCond_noADC(Ag) = Ag_ss_new;
InitCond_noADC(Ag_endo) = Ag_endo_ss_new(end);
InitCond_noADC(Ag_lys) = Ag_lys_ss_new(end);
fprintf("The optimized value for kAgProduction is %e nmol/hr.\n",optimal_AgProduction)

fprintf("Checking final answer... ")
p_test2 = p_test;
p_test2.kAgProduction = optimal_AgProduction;
[T0,Y0] = ode45(eqns_file,[0 time_ss],InitCond_noADC,ode_options,p_test2);

% Plot final figure
if plotAgProdFigure == "Yes"

    % Plot comparison
    figure;
    set(gcf,'color','w')
    yline(Ag_ss,'r--','linewidth',2)
    hold on;
    plot(T0,Y0(:,Ag),'linewidth',2)
    title({'k_{AgProduction} at Steady State = ',[num2str(optimal_AgProduction) ' nmol/hr']},'fontsize',20)
    ylabel('Concentration (nM)','fontsize',20)
    xlabel('Time (hours)','fontsize',20)
    legend(['Ag_s_s',varNames(Ag)])
    set(gca,'FontSize',20)
    grid on; box on;

end

fprintf("Done!\n")

toc