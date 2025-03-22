function [AIC, BIC] = reo_calculate_AICBIC(negative_log_likelihood, params, days)
num_params = length(params);
num_observations = length(days);
log_likelihood_value = -negative_log_likelihood;
[~,~,ic] = aicbic(log_likelihood_value,num_params,num_observations);
AIC = ic.aicc; %corrected AIC
BIC = ic.bic;
end