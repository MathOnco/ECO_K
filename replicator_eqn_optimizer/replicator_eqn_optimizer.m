function finalM = replicator_eqn_optimizer(samples, M, lb, ub)
% replicator_eqn_optimizer: Iteratively optimize the payoff matrix M by removing non-significant interactions.
% This function implements the optimization routine described in the Methods section.
% It minimizes the negative log-likelihood while using AIC/BIC criteria to decide on removing interactions.

% Step 1: Create optimization problem and set up all options
nsp = 80; % number of start points used in optimizer problem
opts = optimoptions(@fmincon, 'Algorithm', 'sqp');
% opts = optimoptions(@fmincon, 'Algorithm', 'interior-point'); % alternative algorithm (commented out)
ms = MultiStart('UseParallel', true);
ms.Display = 'off';

% Initialize records for BIC and AIC during the iterative optimization process
BIC_record = nan(numel(M), 2); % [BIC value, index removed]
AIC_record = nan(numel(M), 2); % [AIC value, index removed]
finalBIC = inf;
finalAIC = inf;
finalM = M;
max_iterations = length(M(M~=0)); % Maximum iterations based on number of non-zero parameters
iteration_count = 0;

% Iteratively remove interactions until no further improvement in BIC is possible
while sum(M(:)~=0) > 2 && iteration_count < max_iterations
    iteration_count = iteration_count + 1;
    M_copy = M;
    idx = find(M_copy ~= 0);
    
    % Optimize the current set of non-zero parameters using the likelihood function
    problem = createOptimProblem('fmincon', 'objective', ...
        @(params) reo_likelihood_function(params, samples, M_copy, []), ...
        'x0', M_copy(idx), 'lb', lb(idx), 'ub', ub(idx), 'options', opts);
    rs = RandomStartPointSet('NumStartPoints', nsp);
    points = list(rs, problem);
    [estimated_params, negative_log_likelihood] = run(ms, problem, CustomStartPointSet(points));
    
    % Calculate AIC and BIC for the current model (see Methods for the log-likelihood model)
    days = unique([samples{2,:}]);
    [AIC_all, BIC_all] = reo_calculate_AICBIC(negative_log_likelihood, estimated_params, days);
    
    % Update final model if current BIC is lower
    if BIC_all < finalBIC
        finalBIC = BIC_all;
        finalAIC = AIC_all;
        finalM = M;
    end

    % Initialize arrays to store AIC and BIC for models with one parameter removed
    BIC_allMinusOne = inf(1, nnz(M)); % For each non-zero parameter, compute BIC after its removal
    AIC_allMinusOne = inf(1, nnz(M));

    valid_indices = find(M ~= 0);
    % Loop over each non-zero parameter to test removal
    for k = 1:length(valid_indices)
        i = valid_indices(k);
        M_copy = M;
        M_copy(i) = 0; % Set one interaction parameter to zero
        idx = find(M_copy ~= 0);

        % Update the optimization problem with the new set of parameters after removal
        problem = createOptimProblem('fmincon', 'objective', ...
            @(params) reo_likelihood_function(params, samples, M_copy, []), ...
            'x0', M_copy(idx), 'lb', lb(idx), 'ub', ub(idx), 'options', opts);
        points = list(rs, problem);
        [estimated_params_minus_one, negative_log_likelihood_minus_one] = run(ms, problem, CustomStartPointSet(points));
        [AIC_allMinusOne(k), BIC_allMinusOne(k)] = reo_calculate_AICBIC(negative_log_likelihood_minus_one, estimated_params_minus_one, days);
    end

    % Identify which parameter removal gives the lowest BIC (and AIC)
    [min_BIC, min_index] = min(BIC_allMinusOne);
    [min_AIC, min_AIC_index] = min(AIC_allMinusOne);

    ia = valid_indices(min_index); % Actual index in M for the lowest BIC when removed
    ix = valid_indices(min_AIC_index); % Actual index in M for the lowest AIC when removed

    % If removal of the parameter improves BIC, update M by permanently removing that parameter.
    if min_BIC < finalBIC
        finalBIC = min_BIC;
        finalAIC = min_AIC;
        M(ia) = 0;
        finalM = M;
    else
        % Even if removal does not lower BIC, remove the parameter (iterative removal continues)
        M(ia) = 0;
    end
    % Record the BIC and AIC for tracking optimization progress
    BIC_record(sum(M(:) == 0), :) = [min_BIC, ia];
    AIC_record(sum(M(:) == 0), :) = [min_AIC, ix];
end

if iteration_count >= max_iterations
    warning('Max iterations reached, exiting loop to prevent infinite run');
end

end
