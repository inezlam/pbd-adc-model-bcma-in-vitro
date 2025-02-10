% Driver File to generate figures (in vitro anti-BCMA PBD ADC model)
clear all; clc;

%% Select figure number to plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FigureNumber = 9; % FigureNumber can be 2, 3, 4, 5, 6, 7, 8, or 9
% Select which subfigure to plot by using "yes" or "no" below as needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setup files
eqns_file = @BCMA_eqns_SingleCell;
setup_file = 'BCMA_model_setup_SingleCell.m';
ode_options = odeset('MaxStep',5e-1, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-2);

switch FigureNumber

    case 2 
        %% Figure 2 - Warhead-Related Figures
        plotFigure2B = "yes"; % Optimization of KA - Crosslink Formation vs. Dose
        plotFigure2C = "yes"; % Validation - Crosslink Formation vs. Time
        plotFigure2D = "yes"; % Optimization of kkill - Cell Survival vs. Dose
        run('plotFigure2.m');

    case 3
        %% Figure 3 - ADC Optimization and Validation Figures
        % Isotype ADC
        plotFigure3A = "yes"; % Optimization - Cell Survival vs. Dose at 96 hrs
        plotFigure3B = "yes"; % Validation - Cell Survival vs. Dose at 72 hrs
        run('plotFigure3AB.m');

        % Standard ADC
        plotFigure3C = "yes"; % Optimization - Cell Survival vs. Dose at 96 hrs
        plotFigure3D = "yes"; % Validation - Cell Survival vs. Dose at Varying sBCMA
        run('plotFigure3CD.m');

    case 4
        %% Figure 4 - Sensitivity Analysis
        % Local Sensitivity
        plotFigure4A = "yes";
        plotFigure4B = "yes";
        plotFigure4C = "yes";
        plotFigure4D = "yes";
        run('plotFigure4_LocalSA.m');
        
        % Global Sensitivity
        plotFigure4E = "yes";
        run('plotFigure4_GlobalSA.m');

    case 5
        %% Figure 5 - DAR Effects
        plotFigure5A = "yes";
        plotFigure5B = "yes";
        plotFigure5C = "yes";
        plotFigure5D = "yes";
        run('plotFigure5_DAR.m');

    case 6
        %% Figure 6 - Warhead Potency Effects
        plotFigure6A = "yes";
        plotFigure6B = "yes";
        run('plotFigure6_Potency.m');
        
        plotFigure6C = "yes";
        plotFigure6D = "yes";
        plotFigure6E = "yes";
        run('plotFigure6_Lipophilicity.m');

    case 7
        %% Figure 7 - Warhead Tracking Plots
        plotFigure7B = "yes";
        plotFigure7C = "yes";
        run('plotFigure7_WarheadTracking.m');

    case 8
        %% Figure 8 - Volume Plots
        plotFigure8A = "yes";
        plotFigure8B = "yes";
        run('plotFigure8_Volume.m');

    case 9
        %% Figure 9 - ADC Design Heatmaps
        plotFigure9 = "yes";
        run('plotFigure9_ADCDesign.m');

end














