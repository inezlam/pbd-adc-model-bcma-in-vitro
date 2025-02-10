function cost_out = costfxn_AgProd2(fxn, options, x0, Ag_ss, p, time, AgSpecies, AgEndoSpecies, AgLysSpecies, InitCond, numTimePoints)
%{ 
This cost function takes in the parameters to be optimized along with the
experimental time points and concentrations. It also takes in the parameter 
vector for the particular patient since that varies.

Because lsqnonlin requires the function to only have the parameters to be 
optimized as input, this cost function is used in an anonymous function in the
driver file to create the input function for lsqnonlin. 

Output is a vector of the difference between the model concentration and the 
experimental concentration at every experimental time point. This output is
minimized by lsqnonlin by adjusting the parameter values in p. 
%}

%% Parameters
p.kAgProduction = x0;
% timePoints = linspace(0,time,numTimePoints);
timePoints = linspace(time/2,time,numTimePoints);

 % Get dose response curve
% [X, Y] = dose_response_v2(fxn, p, time, doses, Species, y0);
[T,Ag_ss_calculated,~,~] = Ag_ss_calculation2(fxn,options,p,time,AgSpecies,AgEndoSpecies,AgLysSpecies,InitCond,x0);

% Initialize output vector
cost_out = zeros(1,numTimePoints);

% For each experimental time point
% cost_out = 
for i=1:numTimePoints
    
    % Finding the differences between the output time and the experimental time
    t_diff = abs(T-timePoints(i));
    
    % Finding the value of T that minimizes the difference with the experimental value 
    [~, t_index] = min(t_diff);
    
    % Comparing the output concentration to the experimental concentration at the closest
    % matching times
    % The goal of the optimization is to minimize this value for each experimental point
    cost_out(i) = Ag_ss_calculated(t_index) - Ag_ss; %Y(x_index) - y_exp(i);
end

end
