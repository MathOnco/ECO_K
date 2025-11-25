function [bootstrap_estimates, num_significant,num_significant2,num_significant3] = bootstrap_func(finalS, samples, num_bootstrap, ub, lb)
    fprintf('Bootstrapping has begun at: %s\n', datetime('now'));
    num_samples_per_bootstrap = round(num_bootstrap / size(samples, 2));
    num_bootstrap = size(samples, 2) * num_samples_per_bootstrap;
    bootstrap_estimates = zeros(num_bootstrap, length(finalS(finalS ~= 0)));
    opts = optimoptions(@fmincon, 'Display', 'off', 'Algorithm', 'sqp');

    % Set up the MultiStart object
    ms = MultiStart('UseParallel', false, 'Display', 'off');

    % Create a cell array to hold the futures
    futures = cell(num_bootstrap, 1);

    % Loop to submit the optimization problems to the parallel pool
    index = 1;
    for p = 1:size(samples, 2)
        for i = 1:num_samples_per_bootstrap
            % Resample the data
            resampled_samples = resample_data(samples, p, finalS);
            resampled_samples = resampled_samples(:, p);

            idx = find(finalS ~= 0);

            % Create the optimization problem
            problem = createOptimProblem('fmincon', 'objective', ...
                @(params) likelihood_function(params, resampled_samples, finalS, []), ...
                'x0', finalS(idx), 'lb', lb(idx), 'ub', ub(idx), 'options', opts);

            % Use parfeval to submit the job to the parallel pool
            futures{index} = parfeval(@runMultiStart, 1, ms, problem);
            index = index + 1;
        end
    end

    % Retrieve the results
    for i = 1:num_bootstrap
        bootstrap_estimates(i, :) = fetchOutputs(futures{i});
    end

    % Analyze bootstrap results
    [param_means, param_std_err, p_values, ci_lower, ci_upper, column_diffs, num_significant, num_significant2, num_significant3] = analyzeBootstrap(bootstrap_estimates, finalS);

    % Plot bootstraps
    plotBootstrapResults(param_means, param_std_err, p_values, ci_lower, ci_upper, bootstrap_estimates, finalS);
end

function [param_means, param_std_err, p_values, ci_lower, ci_upper, column_diffs, num_significant, num_significant2, num_significant3] = analyzeBootstrap(bootstrap_estimates, ~)
    % Calculate means and standard errors of bootstrap estimates
    param_means = mean(bootstrap_estimates);
    param_std_err = std(bootstrap_estimates);

    % Calculate differences between all pairs of columns in bootstrap_estimates
    [m, n] = size(bootstrap_estimates);  % Get the number of rows and columns
    num_pairs = nchoosek(n, 2);  % Number of column pairs
    column_diffs = zeros(m, num_pairs);  % Initialize the column difference matrix

    % Calculate the differences
    counter = 1;
    for i = 1:n
        for j = i+1:n
            column_diffs(:, counter) = bootstrap_estimates(:, i) - bootstrap_estimates(:, j);
            counter = counter + 1;
        end
    end


    % Significance level
    alpha = 0.05;

    % % Number of parameters
    num_params = size(bootstrap_estimates, 2);

    % % Initialize p-values and confidence intervals
    ci_lower = prctile(bootstrap_estimates, (alpha/2)*100);
    ci_upper = prctile(bootstrap_estimates, (1 - alpha/2)*100);

     % Number of parameters
    num_diffs = size(column_diffs, 2);

    % Initialize p-values and confidence intervals
    p_values = zeros(1, num_params);
    ci_lower2 = prctile(column_diffs, (alpha/2)*100);
    ci_upper2 = prctile(column_diffs, (1 - alpha/2)*100);

    % Test for significant difference from zero
    num_significant2 = 0;
    for i = 1:num_diffs
        % If zero is outside the confidence interval, the parameter is significant
        if ci_lower2(i) > 0 || ci_upper2(i) < 0
            num_significant2 = num_significant2 + 1;
        end
    end

    num_significant = 0;
    for i = 1:num_params
         if ci_lower(i) > 0 || ci_upper(i) < 0
            num_significant = num_significant + 1;
        end
        % Compute p-value using a simple t-test approximation
        % Null hypothesis: mean is zero
        t_stat = abs(param_means(i)) / (param_std_err(i) / sqrt(size(bootstrap_estimates, 1)));
        p_values(i) = 2 * (1 - tcdf(t_stat, size(bootstrap_estimates, 1) - 1));  % Two-tailed test
    end

    num_significant3 = sum(p_values<0.05);

end


function plotBootstrapResults(~, ~, ~, ci_lower, ci_upper, bootstrap_estimates, S)

% Plotting the histograms with confidence intervals
num_params = length(S(S~=0));
num_cols = ceil(num_params / 3); % Calculate the number of columns needed

num_bootstraps=length(bootstrap_estimates);
figure();
for i = 1:num_params
    subplot(3, num_cols, i); % Create a subplot for each parameter
    histogram(bootstrap_estimates(:, i), ceil(2 * num_bootstraps^(1/3)));
    hold on;
    y_limits = get(gca, 'YLim');
    plot([ci_lower(i) ci_lower(i)], y_limits, 'r--', 'LineWidth', 2); % lower CI
    plot([ci_upper(i) ci_upper(i)], y_limits, 'r--', 'LineWidth', 2); % upper CI
    hold off;
    title(sprintf('Parameter %d', i));
    xlabel('Parameter Value');
    ylabel('Frequency');
end
end

% Resample data for a specific sample
function resampled_samples = resample_data(samples, sample_idx, payoff_matrix)
resampled_samples = cell(size(samples,1),1);
days = samples{2, sample_idx};
freqs = samples{1, sample_idx};
[~, Y] = ode45(@(t, y) replicatorEqn(t, y, payoff_matrix), days, mean(freqs(:, 1:3), 2)); % this term 'mean(freqs(:, 1:3), 2)' is the guess of initial clone sizes, mean of first 3 observations
num_samples = length(samples{2, sample_idx});
num_samples = max(num_samples, 4);
resample_indices = sort(randperm(length(days), num_samples-1));
% resampled_freqs = freqs(:, resample_indices);
Y=Y';
resampled_freqs = Y(:, resample_indices);
resampled_days = days(resample_indices);
resampled_samples{1, sample_idx} = resampled_freqs;
resampled_samples{2, sample_idx} = resampled_days;
end

function est_params = runMultiStart(ms, problem)
    [est_params, ~] = run(ms, problem, 80);
end