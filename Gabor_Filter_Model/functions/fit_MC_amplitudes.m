function [parameters, fval] = fit_MC_amplitudes(activation_fun, gaborResponses, dataToFitModel, max_iter, use_SurrSupp, raw_gabors, fit_outerSurroundLimit)
%% Now the modeling part:
% here I define how to call my function to set the dataToFit and gaborResponses:
fun = @(x) model_function(x, activation_fun, dataToFitModel, gaborResponses, use_SurrSupp);

% Fix for reproducibility:
rng("default");

if ~exist("max_iter", "var"); max_iter = 100; end
if ~exist("raw_gabors", "var"); raw_gabors = 0; end

maxSurroundDistance = ceil((sqrt(power(size(gaborResponses, 1), 2)+power(size(gaborResponses, 2), 2)) ./ 2) ./ (160./130));

all_fvals = nan(max_iter, 1);
for iter = 1:max_iter
    maxEvals = 10000;
    % options = optimset("Display", "none", "MaxFunEvals", maxEvals);
    % options = optimset("Display", "off", "MaxFunEvals", maxEvals, "TolX", 1e-10, "TolFun", 1e-10);
    options = optimset("Display", "iter", "MaxFunEvals", maxEvals, "TolX", 1e-20, "TolFun", 1e-20);
    if raw_gabors
        if fit_outerSurroundLimit
            if use_SurrSupp == 1  % This is the new version with divisive surround suppression!
%                 % firstGuess = [15 1 0.1, 21];
%     %             firstGuess = [74 1 0.1, 75];
%                 firstGuess = [maxSurroundDistance 1 0.1, maxSurroundDistance];
%                 firstGuess = rand(size(firstGuess)) .* firstGuess;  % I scale the random guesses to reasonable ranges
%                 % [parameters, fval] = fminsearchbnd(fun, firstGuess, [7.5 0 -inf 7.5], [22 inf inf 22], options);
%                 [parameters, fval] = fminsearchbnd(fun, firstGuess, [0 0 -inf 0], [maxSurroundDistance inf inf maxSurroundDistance], options);

            
                % THIS IS THE RERUN 2024-08-15 prohibiting the corners of
                % the 37x37 sourround grid!
                firstGuess = [14 1 0.1, 15, 1];
                firstGuess = rand(size(firstGuess)) .* firstGuess;  % I scale the random guesses to reasonable ranges

                if 0  % No general_offset but center-weight:
                    [parameters, fval] = fminsearchbnd(fun, firstGuess, [7.5 0 0 8.5 0], [14 1 0 15 1], options);
                elseif 0  % general_offset and center-weight:
                    [parameters, fval] = fminsearchbnd(fun, firstGuess, [7.5 0 -inf 8.5 0], [14 1 inf 15 1], options);
                elseif 0  % neither general_offset nor center-weight:
                    [parameters, fval] = fminsearchbnd(fun, firstGuess, [7.5 0 0 8.5 1], [14 1 0 15 1], options);
                else  % general_offset, but no center-weight (the original version):
                    % %% responses = add_surround_supression(gabor_responses, sourround_threshold, surround_weight, surround_offset, surround_threshold_outer, surroundMode, center_weight)
%                     [parameters, fval] = fminsearchbnd(fun, firstGuess, [7.5 0 -inf 8.5 1], [14 1 inf 15 1], options);
% %                     [parameters, fval] = fminsearchbnd(fun, firstGuess, [7.5 0 0 8.5 1], [14 1 1 15 1], options);
                    [parameters, fval] = fminsearchbnd(fun, firstGuess, [7.5 0 -inf 8.5 0], [14 1 inf 15 inf], options);
                end
            elseif use_SurrSupp == 2  % This is the new version with divisive surround suppression!
                % SurroundSupression linear fit
                firstGuess = [14 1 1 15, 1];  % I scale the random guesses to reasonable ranges
                firstGuess = rand(size(firstGuess)) .* firstGuess;  % I scale the random guesses to reasonable ranges

                if 0  % No general_offset but center-weight:
                    [parameters, fval] = fminsearchbnd(fun, firstGuess, [7.5 0 0 8.5 0], [14 inf 0 15 inf], options);
                elseif 0  % general_offset and center-weight:
                    [parameters, fval] = fminsearchbnd(fun, firstGuess, [7.5 0 -inf 8.5 0], [14 inf inf 15 inf], options);
                elseif 0  % neither general_offset nor center-weight:
                    [parameters, fval] = fminsearchbnd(fun, firstGuess, [7.5 0 0 8.5 1], [14 inf 0 15 1], options);
                else  % general_offset, but no center-weight (the original version):
%                     [parameters, fval] = fminsearchbnd(fun, firstGuess, [7.5 0 -inf 8.5 1], [14 inf inf 15 1], options);
% %                     [parameters, fval] = fminsearchbnd(fun, firstGuess, [7.5 0 0 8.5 1], [14 inf inf 15 1], options);
                    [parameters, fval] = fminsearchbnd(fun, firstGuess, [7.5 0 -inf 8.5 0], [14 inf inf 15 inf], options);
                end
            end
        else
            firstGuess = [15 1 0.1];
            firstGuess = rand(size(firstGuess)) .* firstGuess;  % I scale the random guesses to reasonable ranges
            % [parameters, fval] = fminsearch(fun, firstGuess, options);
            [parameters, fval] = fminsearchbnd(fun, firstGuess, [7.5 0 -inf], [22 inf inf], options);
        end
    else
        if use_SurrSupp == 1
            try
                % SurroundSupression nonlinear fit
                firstGuess = rand(1, 6) .* [1 400 0.1 0.1 15 1];  % I scale the random guesses to reasonable ranges
                % [parameters, fval] = fminsearchbnd(fun, firstGuess, [0 0 -inf -inf 7.5 -1], [inf inf inf inf 22 inf], options);
                % [parameters, fval] = fminsearchbnd(fun, firstGuess, [0 0 -inf 0 7.5 -1], [inf inf inf 0 22 inf], options);
                [parameters, fval] = fminsearchbnd(fun, firstGuess, [0 0 -inf 0 7.5 -1], [inf inf inf 0 15 inf], options);  % Changed the outer limit to 15 for the second round of Review 2024-08-14
            catch
                % SurroundSupression linear fit
                firstGuess = rand(1, 4) .* [1 1 15 1];  % I scale the random guesses to reasonable ranges
                % [parameters, fval] = fminsearchbnd(fun, firstGuess, [0 0 -inf -inf 7.5 -1], [inf inf inf inf 22 inf], options);
                % [parameters, fval] = fminsearchbnd(fun, firstGuess, [-inf -inf 7.5 -1], [inf inf 22 inf], options);
                [parameters, fval] = fminsearchbnd(fun, firstGuess, [-inf -inf 7.5 -1], [inf inf 15 inf], options);  % Changed the outer limit to 15 for the second round of Review 2024-08-14
            end
        elseif use_SurrSupp == 0
            try
                % nonlinear without surround-supression
                % firstGuess = rand(1, 4) .* [1 10 0 0];  % I scale the random guesses to reasonable ranges
                firstGuess = rand(1, 4) .* [1 400 0 0];  % I scale the random guesses to reasonable ranges
                
                % Alternatively use a bounded version:
                % This is a custom function from FileExchange, by: John D'Errico 2006.
                % [parameters, fval] = fminsearchbnd(fun, firstGuess, [0 0 -inf -inf], [inf inf inf inf], options);
                [parameters, fval] = fminsearchbnd(fun, firstGuess, [0 0 -inf 0], [inf inf inf 0], options);
            %     [parameters, fval] = fminsearch(fun, firstGuess, options);
            catch ME
                % Linear without surround-supression
                disp(ME)
                firstGuess = rand(1, 2) .* [1 1];  % I scale the random guesses to reasonable ranges
                [parameters, fval] = fminsearchbnd(fun, firstGuess, [-inf -inf], [inf inf], options);
            end
        end
    end


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


