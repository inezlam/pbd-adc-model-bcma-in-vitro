function survival_out = GS_costfxn(fxn, x_exp, y_exp, p, time, doses, Species, y0)
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

% Translate inputs
%% Parameters
% 
% p.ki = x0(1);
% p.frec = x0(2);
% p.ke = x0(3);
% p.kiAg = x0(4);
% p.frecAg = x0(5);
% p.keAg = x0(6);

 % Get dose response curve
[X, Y] = dose_response(fxn, p, time, doses, Species, y0);

% Initialize output vector
survival_out = zeros(1,length(x_exp));

% For each experimental time point
for i=1:length(x_exp)
    
    % Finding the differences between the output time and the experimental time
    x_diff = abs(X-x_exp(i));
    
    % Finding the value of T that minimizes the difference with the experimental value 
    [~, x_index] = min(x_diff);
    
    % Comparing the output concentration to the experimental concentration at the closest
    % matching times
    % The goal of the optimization is to minimize this value for each experimental point
    survival_out(i) = Y(x_index) - y_exp(i);
end

end
