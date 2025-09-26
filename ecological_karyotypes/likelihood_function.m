% Combined likelihood function
function negative_log_likelihood = likelihood_function(params, samples, M, solMat)
negative_log_likelihood = 0;
OUTDIR = 'Results/';
for p = 1:size(samples, 2)
    days = samples{2, p};
    X = samples{1, p};  % observed clone frequencies (X)
    M_matrix = M;
    idx = find(M_matrix ~= 0);
    M_matrix(idx) = params;
    origin_ID = samples{3,1};
    if ~isempty(solMat)
        % Remove 'Origin-' prefix from origin_ID if provided
        if ~isempty(origin_ID)
            originWithoutPrefix = strrep(origin_ID, 'Origin-', '');
            matFileName = ['color_files/clone_colors_', originWithoutPrefix, '.mat'];
        else
            matFileName = '';
        end
        % Attempt to load the color palette
        if ~isempty(matFileName) && exist(matFileName, 'file')
            % Load the cloneColorTable variable from the .mat file
            load(matFileName, 'cloneColorTable');
            % Filter by clones in this Replicate Group
            idx = ismember(cloneColorTable.CloneID, samples{4, 1});
            cloneColorTable = cloneColorTable(idx, :);
            % Extract the RGB color values into nodecols
            color_map = [cloneColorTable.R, cloneColorTable.G, cloneColorTable.B];
        else
            % Use a default colormap if the .mat file is not found
            disp('Color table file not found. Using MATLAB default colormap.');
            color_map = lines(size(samples{1, 1}, 1));
        end

        x_pred = solMat{1, p};  % predicted clone frequencies (x)
        t = solMat{2, p};       % time vector

        % Calculate the new column such that each row sums to 1
        additional_column = 1 - sum(x_pred, 2);
        referenceClone = readtable(['fitCloneOutput/referenceClones/' samples{5} '_reference_clone.txt'], 'Delimiter', '\t');
        referenceClone = referenceClone.Removed_Clone + 1;
        x_pred = [x_pred(:, 1:referenceClone-1), additional_column, x_pred(:, referenceClone:end)];

        figure();
        set(gca, 'ColorOrder', color_map, 'NextPlot', 'replacechildren');
        plot(t, x_pred, 'LineWidth', 4);
        hold on;

        % Plot actual data
        for i = 1:size(X, 1)
            plot(days, X(i, :), '--o', 'MarkerSize', 12, 'LineWidth', 1, ...
                'Color', color_map(i, :), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', color_map(i, :));
        end

        % Set axis labels and title
        xlabel('Day', 'FontSize', 14, 'FontWeight', 'bold');
        ylabel('Frequency', 'FontSize', 14, 'FontWeight', 'bold');

         title([samples{3} ' Sample ' num2str(samples{6}) ' Replicate ' num2str(p) ' fitClone Model Fit']);

        % Add legend and grid
        legend(cellstr(samples{4, p}), 'Location', 'northeast', 'FontSize', 12);
        grid on;
        set(gca, 'GridColor', [0.9, 0.9, 0.9], 'LineWidth', 0.5);
        hold off;

        g = gcf;
        imageSaveName = [samples{5} '_Replicate_' num2str(p) '_fitCloneModelFit.png'];
        savePlace = ['Results/fitClonePlots/' imageSaveName];
        exportgraphics(g, savePlace, 'Resolution', 300);
        % Find the closest values in t to days and use the indices to index x_pred
        closest_indices = zeros(length(days), 1);
        for i = 1:length(days)
            [~, closest_indices(i)] = min(abs(t - days(i)));
        end
        % Index x_pred using closest_indices
        x_pred = x_pred(closest_indices, :);
        res = X - x_pred';
        err = sum(sum((X - x_pred').^2, 1));
        disp(err)
        % Calculate standard deviation for each compartment
        sigma = sqrt(1/(4*100)); % Using the value N=1146 from the paper
        % Calculate the log-likelihood for each residual
        log_likelihood_matrix = zeros(size(res));
        for j = 1:size(res, 2)
            log_likelihood_matrix(:, j) = log(normpdf(res(:, j), 0, sigma));
        end
        % Sum the log-likelihoods across all compartments and time points
        total_log_likelihood = sum(log_likelihood_matrix(:));
        % Update negative log-likelihood
        negative_log_likelihood = negative_log_likelihood - total_log_likelihood;
    else
        [~, x_pred] = ode45(@(t, y) replicatorEqn(t, y, M_matrix), days, mean(X(:, 1:3), 2)); % initial condition is the mean of the first 3 observations
        [rowsF, colsF] = size(X);
        [rowsY, colsY] = size(x_pred');

        % Check if the number of columns in x_pred' is less than in X
        if colsY < colsF
            % Extract the last row of x_pred
            lastRow = x_pred(end,:);
            % Number of additional rows to add to x_pred so x_pred' will have enough columns
            padRows = colsF - colsY;
            % Create a padding matrix by replicating the last row of x_pred
            padMatrix = repmat(lastRow, padRows, 1);
            % Pad x_pred (adding rows)
            x_pred = [x_pred; padMatrix];
        end
        res = X - x_pred';
        % Calculate standard deviation for each compartment
        sigma = sqrt(1/(4*100)); % Using the value N=1146 from the paper
        % Calculate the log-likelihood for each residual
        log_likelihood_matrix = zeros(size(res));
        for j = 1:size(res, 2)
            log_likelihood_matrix(:, j) = log(normpdf(res(:, j), 0, sigma));
        end
        % Sum the log-likelihoods across all compartments and time points
        total_log_likelihood = sum(log_likelihood_matrix(:));
        % Update negative log-likelihood
        negative_log_likelihood = negative_log_likelihood - total_log_likelihood;
    end
end
end
