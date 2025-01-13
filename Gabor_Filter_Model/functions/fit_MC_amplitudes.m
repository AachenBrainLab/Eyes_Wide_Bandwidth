function [parameters, fval] = fit_MC_amplitudes(activation_fun, gaborResponses, dataToFitModel, max_iter, use_SurrSupp)
%% Now the modeling part:
% here I define how to call my function to set the dataToFit and gaborResponses:
fun = @(x) model_function(x, activation_fun, dataToFitModel, gaborResponses, use_SurrSupp);

% Fix for reproducibility:
rng("default");

if ~exist("max_iter", "var"); max_iter = 100; end

maxSurroundDistance = ceil((sqrt(power(size(gaborResponses, 1), 2)+power(size(gaborResponses, 2), 2)) ./ 2) ./ (160./130));

all_fvals = nan(max_iter, 1);
for iter = 1:max_iter
    maxEvals = 10000;
    % options = optimset("Display", "none", "MaxFunEvals", maxEvals);
    % options = optimset("Display", "off", "MaxFunEvals", maxEvals, "TolX", 1e-10, "TolFun", 1e-10);
    options = optimset("Display", "iter", "MaxFunEvals", maxEvals, "TolX", 1e-20, "TolFun", 1e-20);

    % SurroundSupression linear fit
    firstGuess = [14 1 1 15, 1];  % I scale the random guesses to reasonable ranges
    firstGuess = rand(size(firstGuess)) .* firstGuess;  % I scale the random guesses to reasonable ranges
    [parameters, fval] = fminsearchbnd(fun, firstGuess, [7.5 0 -inf 8.5 0], [14 inf inf 15 inf], options);
    

    if iter == 1
        best_params = parameters;
    end
    if fval < min(all_fvals(1:iter-1))
        best_params = parameters;
    end
    all_fvals(iter) = fval;

    disp(iter + ": " + best_params);
    if mod(iter, 50) == 0
        figure(101)
        plot(sort(all_fvals)); drawnow;
    end
end

parameters = best_params;
fval = min(all_fvals(1:iter-1));



% define cost function to minimize
function prediction_error = model_function(x, activation_fun, dataToFitModel, gabor_responses, use_SurrSupp)
    % Apply the activation function:
    responses = activation_fun(gabor_responses, x, use_SurrSupp);

    % replace possible NaNs (e.g.: div. 0)
    responses(isnan(responses)) = 0;

    % function to minimize the error between observation and prediction
    prediction_error = rms(reshape(dataToFitModel - responses, [], 1));
end


end


