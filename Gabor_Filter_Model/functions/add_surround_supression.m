function responses = add_surround_supression(gabor_responses, surround_threshold, surround_weight, surround_offset, surround_threshold_outer, surroundMode, center_weight)
if ~exist("surround_weight", "var"); surround_weight = 1; end
if ~exist("surround_offset", "var"); surround_offset = 0; end
if ~exist("surround_threshold_outer", "var"); surround_threshold_outer = inf; end
if ~exist("surroundMode", "var"); surroundMode = 1; end
if ~exist("center_weight", "var"); center_weight = 1; end

if surround_threshold < 7.5; surround_threshold = 7.5; end
if surround_threshold > 14; surround_threshold = 14; end
if surround_threshold_outer < surround_threshold + 1; surround_threshold_outer = surround_threshold + 1; end
if surround_threshold_outer > 15; surround_threshold_outer = 15; end

% surround_threshold: in the range of [7.5 sqrt(2*15**2)~21]
% relu_threshold: in principle: [-1 1]  % Most likely: [~0 0.3]
% center_threshold = 15./2;

% set basic parameters for how to computed the Gabor-responses:
centerGaborIDs = floor(size(gabor_responses, 1:2) ./ 2)+1;  % center_gabor_id = 19;
px_per_deg = 160./130;
% dist_mtx = sqrt(power(-18:18, 2) + power(-18:18, 2)') ./ px_per_deg;
dist_mtx = sqrt(power((1:size(gabor_responses, 1))-centerGaborIDs(1), 2)' + power((1:size(gabor_responses, 2))-centerGaborIDs(2), 2)) ./ px_per_deg;
surround_mask = single((dist_mtx >= surround_threshold) & (dist_mtx <= surround_threshold_outer));
surround_mask(surround_mask == 0) = nan;
surround_mask = reshape(surround_mask, size(surround_mask, 1), size(surround_mask, 2), 1, 1, 1, 1);


% center_response = gabor_responses(center_gabor_id, center_gabor_id, SF_ids, :, :, :);
center_response = gabor_responses(centerGaborIDs(1), centerGaborIDs(2), :, :, :, :);
surround_response = nanmean(nanmean(surround_mask .* gabor_responses, 1), 2);

if surroundMode == 1
    responses = (center_weight .* center_response) - (surround_weight .* surround_response) + surround_offset;
elseif surroundMode == 2
    % responses = center_response ./ (1 + surround_weight .* surround_response) + surround_offset;
    % responses = (center_weight .* center_response) ./ (1 + surround_weight .* surround_response) + surround_offset;
    divisor = surround_weight .* (surround_response + surround_offset);
    divisor(divisor < 0) = 0;
    %  divisor(divisor > 1) = 1;  % The divisive suppresion is weird when this is outside of the range [0 1]!
    responses = (center_weight .* center_response) ./ (1 + divisor);
else
    error("By default 'surroundMode' should be one. Otherwise further case are not implemented!")
end
% responses = nanmean(nanmean(center_response - surround_response, 6), 3);

% This step was previously outside of this function. Having this in here
% now replaces the apply SurroundSupression function!
responses(responses < 0) = 0;

% responses = responses + surround_offset




