function [hatM_all, num_significant, num_significant2, num_significant3] = reo_bootstrap(finalS, X, beta, ub, lb)
    % REO_BOOTSTRAP performs the bootstrap routine described in Section "Bootstrapping" of the Methods.
    % 
    %   finalS      - Final parameter solution to the replicator equation (defines payoff matrix M).
    %   X           - A 2 x N cell array (original data), analogous to \mathbf{X} in Ω.
    %   beta        - Number of bootstrap samples (was num_bootstrap).
    %   ub, lb      - Upper and lower bounds for optimization.
    %
    %   hatM_all          - Stores all bootstrap parameter estimates (\hat{M}_1, ..., \hat{M}_β).
    %   num_significant   - Count of parameters whose 95% CI excludes zero.
    %   num_significant2  - Count of column pair-differences excluding zero.
    %   num_significant3  - Count of parameters with p-values < 0.05.

    fprintf('Bootstrapping has begun at: %s\n', datetime('now'));

    % Determine how many bootstrap samples to draw from each column in X.
    num_samples_per_bootstrap = round(beta / size(X, 2));
    beta = size(X, 2) * num_samples_per_bootstrap;  % Re-calc actual beta if rounding changed it

    % Preallocate the array that will hold all bootstrap estimates of the non-zero parameters
    hatM_all = zeros(beta, length(finalS(finalS ~= 0)));

    % Set up the optimization options and MultiStart
    opts = optimoptions(@fmincon, 'Display', 'off', 'Algorithm', 'sqp');
    ms   = MultiStart('UseParallel', false, 'Display', 'off');

    % We'll store asynchronous "futures" here
    futures = cell(beta, 1);

    % Index to track how many total bootstrap problems are queued
    index = 1;

    %====================================================================
    % MAIN BOOTSTRAP LOOP: for each column in X and for each bootstrap
    % sample, we resample and optimize to get \hat{M}_b.
    %====================================================================
    for p = 1:size(X, 2)
        for i = 1:num_samples_per_bootstrap

            %------------------------------------------
            % 1) Generate a bootstrap sample X^*_b
            %    (resample with replacement)
            %------------------------------------------
            resampled_samples = generateBootstrapSample(X, p, finalS);

            % The original code then takes: resampled_samples = resampled_samples(:, p);
            % This is how it slices the 2xN cell array. We keep it exactly for consistency.
            resampled_samples = resampled_samples(:, p);

            % Identify indices of nonzero parameters in finalS
            idx = find(finalS ~= 0);

            %------------------------------------------
            % 2) Create the optimization problem
            %------------------------------------------
            problem = createOptimProblem('fmincon', ...
                'objective', @(params) reo_likelihood_function(params, resampled_samples, finalS, []), ...
                'x0',        finalS(idx), ...
                'lb',        lb(idx), ...
                'ub',        ub(idx), ...
                'options',   opts);

            % Queue the optimization with parfeval
            futures{index} = parfeval(@runMultiStart, 1, ms, problem);
            index = index + 1;
        end
    end

    %====================================================================
    % Retrieve all optimization results for \hat{M}_b
    %====================================================================
    for i = 1:beta
        hatM_all(i, :) = fetchOutputs(futures{i});
    end

    %====================================================================
    % 3) Analyze the bootstrap results to obtain means, std. errors, p-values, etc.
    %====================================================================
    [param_means, param_std_err, p_values, ci_lower, ci_upper, column_diffs, ...
        num_significant, num_significant2, num_significant3] = analyzeBootstrap(hatM_all, finalS);

    %====================================================================
    % 4) (Optional) Visualize the bootstrap distribution via histograms
    %====================================================================
    plotBootstrapResults(param_means, param_std_err, p_values, ci_lower, ci_upper, hatM_all, finalS);
end


%==========================================================================
% ANALYZEBOOTSTRAP: Computes means, stdev, confidence intervals, p-values,
% and significance counts from all bootstrap parameter estimates.
%==========================================================================
function [param_means, param_std_err, p_values, ci_lower, ci_upper, column_diffs, ...
    num_significant, num_significant2, num_significant3] = analyzeBootstrap(bootstrap_estimates, ~)

    % 1) Compute the mean and standard deviation across all bootstrap estimates
    param_means = mean(bootstrap_estimates);
    param_std_err = std(bootstrap_estimates);

    % 2) Prepare to calculate pairwise differences among columns (parameters)
    [m, n] = size(bootstrap_estimates); 
    num_pairs = nchoosek(n, 2); 
    column_diffs = zeros(m, num_pairs); 

    % 3) Fill in the matrix of differences: (col_i - col_j) for all pairs
    counter = 1;
    for i = 1:n
        for j = i+1:n
            column_diffs(:, counter) = bootstrap_estimates(:, i) - bootstrap_estimates(:, j);
            counter = counter + 1;
        end
    end

    % 4) Define alpha for significance
    alpha = 0.05;

    % 5) Calculate percentile-based 95% CIs for each parameter
    ci_lower = prctile(bootstrap_estimates, (alpha/2)*100);
    ci_upper = prctile(bootstrap_estimates, (1-alpha/2)*100);

    % 6) Similarly, percentile-based CIs for the pairwise differences
    num_diffs = size(column_diffs, 2);
    ci_lower2 = prctile(column_diffs, (alpha/2)*100);
    ci_upper2 = prctile(column_diffs, (1 - alpha/2)*100);

    %---------------------------------------------------
    % Count how many parameter differences exclude zero
    %---------------------------------------------------
    num_significant2 = 0;
    for i = 1:num_diffs
        if ci_lower2(i) > 0 || ci_upper2(i) < 0
            num_significant2 = num_significant2 + 1;
        end
    end

    %---------------------------------------------------
    % Count how many parameters have CIs excluding zero
    %---------------------------------------------------
    num_significant = 0;
    p_values = zeros(1, n);
    for i = 1:n
        if ci_lower(i) > 0 || ci_upper(i) < 0
            num_significant = num_significant + 1;
        end

        %---------------------------------------------------
        % Compute p-value vs 0 using a simple t-test approach
        %---------------------------------------------------
        t_stat = abs(param_means(i)) / (param_std_err(i) / sqrt(size(bootstrap_estimates, 1)));
        p_values(i) = 2 * (1 - tcdf(t_stat, size(bootstrap_estimates, 1) - 1)); 
    end

    %---------------------------------------------------
    % Count how many parameters have p < 0.05
    %---------------------------------------------------
    num_significant3 = sum(p_values < 0.05);
end


%==========================================================================
% PLOTBOOTSTRAPRESULTS: Plots histograms of the bootstrap distribution for
% each parameter and overlays the 95% CI lines.
%==========================================================================
function plotBootstrapResults(~, ~, ~, ci_lower, ci_upper, bootstrap_estimates, S)

    % Count how many actual parameters are non-zero
    num_params = length(S(S ~= 0));
    num_cols = ceil(num_params / 3); 
    num_bootstraps = size(bootstrap_estimates,1);

    figure();
    for i = 1:num_params
        subplot(3, num_cols, i);
        histogram(bootstrap_estimates(:, i), ceil(2 * num_bootstraps^(1/3)));
        hold on;
        y_limits = get(gca, 'YLim');
        plot([ci_lower(i) ci_lower(i)], y_limits, 'r--', 'LineWidth', 2); % lower bound
        plot([ci_upper(i) ci_upper(i)], y_limits, 'r--', 'LineWidth', 2); % upper bound
        hold off;
        title(sprintf('Parameter %d', i));
        xlabel('Value');
        ylabel('Frequency');
    end
end


%==========================================================================
% GENERATEBOOTSTRAPSAMPLE: Resamples (with replacement) from the solution
% to the replicator equation. Returns a cell array with frequencies & days.
%==========================================================================
function resampled_samples = generateBootstrapSample(X, sample_idx, payoff_matrix)
    % We keep the logic of the original "resample_data" function exactly:
    %   1) Evaluate replicatorEqn with payoff_matrix to get Y.
    %   2) Resample from Y's time points with replacement.
    %
    % Return a 2xN cell array (like the original code does).
    
    resampled_samples = cell(size(X,1), 1);  % same shape as original
    days = X{2, sample_idx};
    freqs = X{1, sample_idx};

    % Solve replicatorEqn from an initial guess (mean of first 3 freq observations)
    [~, Y] = ode45(@(t, y) replicatorEqn(t, y, payoff_matrix), days, mean(freqs(:, 1:3), 2));
    Y = Y';  % to keep dimension consistency

    % Number of time points to sample
    num_samples = length(days);
    num_samples = max(num_samples, 4);  % same logic as the original code

    % Randomly sample indices with replacement
    resample_indices = sort(randperm(length(days), num_samples-1));

    % Collect frequency/time from those indices
    resampled_freqs = Y(:, resample_indices);
    resampled_days  = days(resample_indices);

    % Assign into the cell structure
    resampled_samples{1, sample_idx} = resampled_freqs;
    resampled_samples{2, sample_idx} = resampled_days;
end


%==========================================================================
% RUNMULTISTART: Helper to run the MultiStart optimization (80 starts).
%==========================================================================
function est_params = runMultiStart(ms, problem)
    [est_params, ~] = run(ms, problem, 80);
end
