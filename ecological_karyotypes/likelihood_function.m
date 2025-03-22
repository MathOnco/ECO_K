function negative_log_likelihood = reo_likelihood_function(params, samples, M, solMat)
% reo_likelihood_function: Computes the negative log-likelihood for the replicator equation model.
% This function compares observed clone frequencies with those predicted by the replicator equation,
% incorporating the fitness interactions specified by the payoff matrix M.
% It corresponds to the likelihood model described in the Methods section.

negative_log_likelihood = 0;
OUTDIR = '../../Results/shahDataFits/';

% Loop over each replicate sample
for p = 1:size(samples, 2)
    days = samples{2, p};    % Time points (days)
    X = samples{1, p};       % Observed clone frequencies (data)
    M_matrix = M;            % Current payoff matrix
    idx = find(M_matrix ~= 0);
    M_matrix(idx) = params;  % Insert optimized parameters into the payoff matrix
    
    origin_ID = samples{3,1};
    
    % If a precomputed solution (solMat) is provided, use it for plotting and likelihood calculation
    if ~isempty(solMat)
        % Attempt to load a color palette for plotting based on the origin ID
        if ~isempty(origin_ID)
            originWithoutPrefix = strrep(origin_ID, 'Origin-', '');
            matFileName = [OUTDIR, 'clone_colors_', originWithoutPrefix, '.mat'];
        else
            matFileName = '';
        end
        % Load color palette if available, otherwise use a default colormap
        if ~isempty(matFileName) && exist(matFileName, 'file')
            load(matFileName, 'cloneColorTable');
            idx = ismember(cloneColorTable.CloneID, samples{4, 1});
            cloneColorTable = cloneColorTable(idx, :);
            color_map = [cloneColorTable.R, cloneColorTable.G, cloneColorTable.B];
        else
            disp('Color table file not found. Using MATLAB default colormap.');
            color_map = lines(size(samples{1, 1}, 1));
        end

        % Retrieve precomputed solution: predicted clone frequencies and corresponding time vector
        x_pred = solMat{1, p};  % predicted clone frequencies
        t = solMat{2, p};       % time vector

        % Adjust predicted frequencies so that each row sums to 1 by adding an additional column
        additional_column = 1 - sum(x_pred, 2);
        referenceClone = readtable(['../../Data/shahData/referenceClones/' samples{5} '_reference_clone.txt'], 'Delimiter', '\t');
        referenceClone = referenceClone.Removed_Clone + 1;
        x_pred = [x_pred(:, 1:referenceClone-1), additional_column, x_pred(:, referenceClone:end)];

        % Plot the predicted clone dynamics along with the observed data for visual comparison
        figure();
        set(gca, 'ColorOrder', color_map, 'NextPlot', 'replacechildren');
        plot(t, x_pred, 'LineWidth', 4);
        hold on;
        for i = 1:size(X, 1)
            plot(days, X(i, :), '--o', 'MarkerSize', 12, 'LineWidth', 1, ...
                'Color', color_map(i, :), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', color_map(i, :));
        end
        xlabel('Day', 'FontSize', 14, 'FontWeight', 'bold');
        ylabel('Frequency', 'FontSize', 14, 'FontWeight', 'bold');
        title([samples{3} ' Sample ' num2str(samples{6}) ' Replicate ' num2str(p) ' fitClone Model Fit']);
        legend(cellstr(samples{4, p}), 'Location', 'northeast', 'FontSize', 12);
        grid on;
        set(gca, 'GridColor', [0.9, 0.9, 0.9], 'LineWidth', 0.5);
        hold off;
        g = gcf;
        imageSaveName = [samples{5} '_Replicate_' num2str(p) '_fitCloneModelFit.png'];
        savePlace = ['../../Results/shahDataFits/' imageSaveName];
        exportgraphics(g, savePlace, 'Resolution', 300);
        
        % Match predicted frequencies (x_pred) with the observed time points (days)
        closest_indices = zeros(length(days), 1);
        for i = 1:length(days)
            [~, closest_indices(i)] = min(abs(t - days(i)));
        end
        x_pred = x_pred(closest_indices, :);
        
        % Calculate residuals between observed and predicted clone frequencies
        res = X - x_pred';
        err = sum(sum((X - x_pred').^2, 1));
        disp(err)
        
        % Compute standard deviation for residuals (adding eps to avoid division by zero)
        sigma = std(res, 0, 2) + eps;
        % Calculate log-likelihood for each residual assuming a normal distribution (as per Methods)
        log_likelihood_matrix = zeros(size(res));
        for j = 1:size(res, 2)
            log_likelihood_matrix(:, j) = log(normpdf(res(:, j), 0, sigma));
        end
        % Sum the log-likelihoods across compartments and time points
        total_log_likelihood = sum(log_likelihood_matrix(:));
        negative_log_likelihood = negative_log_likelihood - total_log_likelihood;
    else
        % If no precomputed solution is provided, solve the replicator equation ODE using ode45.
        % The replicator equation is defined in the Methods (see Eqn. for replicator dynamics and the replicatorEqn function)
        [~, x_pred] = ode45(@(t, y) replicatorEqn(t, y, M_matrix), days, mean(X(:, 1:3), 2)); % initial condition is the mean of the first 3 observations
        [rowsF, colsF] = size(X);
        [rowsY, colsY] = size(x_pred');

        % If the ODE solution has fewer time points than observed data, pad the solution accordingly.
        if colsY < colsF
            lastRow = x_pred(end,:);
            padRows = colsF - colsY;
            padMatrix = repmat(lastRow, padRows, 1);
            x_pred = [x_pred; padMatrix];
        end
        res = X - x_pred';
        sigma = std(res, 0, 2) + eps;
        log_likelihood_matrix = zeros(size(res));
        for j = 1:size(res, 2)
            log_likelihood_matrix(:, j) = log(normpdf(res(:, j), 0, sigma));
        end
        total_log_likelihood = sum(log_likelihood_matrix(:));
        negative_log_likelihood = negative_log_likelihood - total_log_likelihood;
    end
end
end

