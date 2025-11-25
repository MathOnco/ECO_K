% REO_CALCULATE_AICBIC
function [AIC, BIC] = calculate_AICBIC(negative_log_likelihood, params, days)
% Returns AIC (actually AICc, corrected AIC) and BIC without using aicbic.
%
% Inputs
%   negative_log_likelihood : scalar NLL at the MLE (−log L)
%   params                  : vector of fitted parameters (used for k = numel(params))
%   days                    : vector of observation indices/times (used for n = numel(days))
%
% Outputs
%   AIC : AICc (corrected AIC)
%   BIC : Bayesian Information Criterion
    k = numel(params);
    n = numel(days);

    % if ~isscalar(negative_log_likelihood) || ~isfinite(negative_log_likelihood)
    %     error('negative_log_likelihood must be a finite scalar.');
    % end
    % if n <= 0 || k < 0
    %     error('Invalid sizes: n must be > 0 and k must be >= 0.');
    % end

    logL = -negative_log_likelihood;

    % Standard AIC
    AIC_std = 2*k - 2*logL;

    % Corrected AIC (AICc). If n - k - 1 <= 0, AICc is undefined; return NaN.
    if n - k - 1 > 0
        AIC = AIC_std + (2*k*(k + 1)) / (n - k - 1);
    else
        % fprintf('Unable to calculate AICc, returning standard AIC.\n');
        AIC = AIC_std;  % or set to Inf if you prefer to always penalize
    end

    % BIC (a.k.a. Schwarz criterion)
    BIC = log(n)*k - 2*logL;
end