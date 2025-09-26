function [finalM, finalBIC] = ecological_karyotypes(samples, M_initial, lb_M, ub_M, beam_width)
% Performs model selection using a Beam Search algorithm to make it less greedy.
% Instead of following the single best path, it explores the 'beam_width' best paths at each step.

%% --- Step 1: Initialization ---
fprintf('Initializing Beam Search with a beam width of %d...\n', beam_width);

% Optimization settings (reused from previous script)
nsp = 80 ;
opts = optimoptions(@fmincon, 'Algorithm', 'sqp', 'Display', 'off');
ms = MultiStart('UseParallel', true, 'Display', 'off');

% --- Initialize the Overall Best Model Tracker ---
% This struct will hold the single best model found across all paths and iterations.
best_overall_model.M = M_initial;
best_overall_model.BIC = inf;

% --- Create the Initial Beam ---
% The search starts with a single "beam" containing the full, initial model.
% A beam is a struct array where each element represents a candidate model.
initial_model.M = M_initial;
[~, initial_model.BIC, ~] = optimize_and_get_bic(initial_model, samples, lb_M, ub_M, ms, opts, nsp);

current_beam = initial_model;
best_overall_model = initial_model; % The initial model is the best one so far

fprintf('Initial full model BIC: %f\n', initial_model.BIC);

% --- Loop Control ---
max_iterations = nnz(M_initial);
min_params = 3; % Stop when models have 3 or fewer parameters
iteration_count = 0;

%% --- Step 2: Main Beam Search Loop ---
while iteration_count < max_iterations
    iteration_count = iteration_count + 1;
    fprintf('\n================== Iteration %d ==================\n', iteration_count);
    
    % --- Expansion Step: Generate all possible child models from the current beam ---
    child_candidates = {}; % Use a cell array to collect all children
    
    for i = 1:length(current_beam)
        parent_model = current_beam(i);
        
        % Find all valid parameters to remove from this parent
        valid_indices_m = find(parent_model.M ~= 0);
        
        % Stop expanding a path if it has too few parameters
        if length(valid_indices_m) <= min_params
            continue;
        end
        
        % Test removing each M parameter
        for k = 1:length(valid_indices_m)
            child_model = parent_model;
            child_model.M(valid_indices_m(k)) = 0;
            [~, child_model.BIC, ~] = optimize_and_get_bic(child_model, samples, lb_M, ub_M, ms, opts, nsp);
            

            % Store how this child was created  <-- ADD THIS LINE
            child_model.description = sprintf('Parent %d: Removed M at index %d', i, valid_indices_m(k));
            
            child_candidates{end+1} = child_model;
        end
        
    end
    
    if isempty(child_candidates)
        fprintf('No further simplifications possible. Terminating search.\n');
        break;
    end
    
    % --- Pruning Step: Select the top 'k' candidates to form the next beam ---
    % Convert cell array to a struct array for easy sorting
    candidate_structs = [child_candidates{:}];
    
    % Get all unique models to avoid redundant paths in the beam
    [~, unique_indices] = unique(arrayfun(@(s) gencode(s.M), candidate_structs, 'UniformOutput', false));
    unique_candidates = candidate_structs(unique_indices);
    
    % Sort all unique candidates by their BIC score
    [~, sort_order] = sort([unique_candidates.BIC]);
    sorted_candidates = unique_candidates(sort_order);
    
    % The new beam is the top 'beam_width' models from the sorted list
    num_to_keep = min(beam_width, length(sorted_candidates));
    next_beam = sorted_candidates(1:num_to_keep);
    
    % --- Update the overall best model tracker FIRST ---
    if ~isempty(next_beam) && (next_beam(1).BIC < best_overall_model.BIC)
        best_overall_model = next_beam(1);
    end

    % --- MODIFIED: More detailed printout at the end of each iteration ---
    % fprintf('\n================== Iteration %d Summary ==================\n', iteration_count);
    % fprintf('Parent models (start of iteration):\n');
    % for p_idx = 1:length(current_beam)
    %     fprintf('  Parent %d -> BIC: %f\n', p_idx, current_beam(p_idx).BIC);
    % end
    % fprintf('\n--- BIC Test Results for This Iteration (%d unique children generated) ---\n', length(sorted_candidates));
    % % Loop through and display the results for all generated models
    % for k = 1:length(sorted_candidates)
    %     candidate = sorted_candidates(k);
    %     fprintf('Test: %-35s -> BIC: %f', candidate.description, candidate.BIC);
    %     if k <= num_to_keep
    %         fprintf('  <-- SELECTED for next beam\n');
    %     else
    %         fprintf('\n');
    %     end
    % end
    % fprintf('--------------------------------------------------------------------\n\n');
    % fprintf('State for NEXT iteration (new beam with %d model(s)):\n', length(next_beam));
    % for b_idx = 1:length(next_beam)
    %     fprintf('--- Beam Model %d (BIC: %f) ---\n', b_idx, next_beam(b_idx).BIC);
    %     disp(next_beam(b_idx).M);
    % end
    % fprintf('\nBest model found so far (lowest BIC):\n');
    % disp('finalM matrix:');
    % disp(best_overall_model.M);
    % fprintf('Lowest BIC: %f\n', best_overall_model.BIC);
    % fprintf('==================================================================\n\n');
    
    current_beam = next_beam;
    
    % fprintf('Generated %d unique child models. New beam has %d models.\n', length(unique_candidates), length(current_beam));
    % fprintf('Best BIC in current beam: %f\n', current_beam(1).BIC);
    % fprintf('Overall best BIC found so far: %f\n', best_overall_model.BIC);

    1+1;

end

%% --- Step 3: Finalization ---
fprintf('\n================== Search Complete ==================\n');
fprintf('The best model found has a BIC of: %f\n', best_overall_model.BIC);

finalM = best_overall_model.M;
finalBIC = best_overall_model.BIC;

end


% --- Helper function for optimization ---
function [nll, bic, estimated_params] = optimize_and_get_bic(model, samples, lb_M, ub_M, ms, opts, nsp)
    % This helper runs the optimization and BIC calculation for a given model structure.
    M_copy = model.M;
    
    idx_m = find(M_copy ~= 0);
    x0 = M_copy(idx_m);
    
    
    lb = lb_M(idx_m);
    ub = ub_M(idx_m);
    
    problem = createOptimProblem('fmincon', 'objective', ...
        @(params) likelihood_function(params, samples, M_copy, []), ...
        'x0', x0, 'lb', lb, 'ub', ub, 'options', opts);
    
    rs = RandomStartPointSet('NumStartPoints', nsp);
    points = list(rs, problem);
    [estimated_params, nll] = run(ms, problem, CustomStartPointSet(points));
    
    days = unique([samples{2,:}]);
    [~, bic] = calculate_AICBIC(nll, estimated_params, days);
end

% --- Helper function to generate a unique code for a model ---
function code = gencode(M)
    % Creates a unique string representation of a model's structure for finding unique models.
    M_binary = M(:) ~= 0;
    code = num2str(M_binary');
end