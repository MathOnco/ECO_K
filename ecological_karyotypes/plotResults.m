function [error] = reo_plotResults(payoff_matrix, samples, OUTDIR, origin_ID, sample_ID, solMat)
    % Initialize error
    error = nan(size(samples, 2), 1);

    % Check if optional parameters are provided
    if nargin < 3 || isempty(OUTDIR)
        OUTDIR = ''; % Do not save images if OUTDIR is empty
    end
    if nargin < 4
        origin_ID = '';
    end
    if nargin < 5
        sample_ID = '';
    end

    % Remove 'Origin-' prefix from origin_ID if provided
    if ~isempty(origin_ID)
        originWithoutPrefix = strrep(origin_ID, 'Origin-', '');
        matFileName = [OUTDIR, 'clone_colors_', originWithoutPrefix, '.mat'];
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

    % Handle solMat or simulate data with payoff_matrix
    if ~isempty(solMat)
        for p = 1:size(samples, 2)
            Y = solMat{1, p};
            T = solMat{2, p};
            freqs = samples{1, p};
            days = samples{2, p};

            % Calculate additional column to ensure each row sums to 1
            additional_column = 1 - sum(Y, 2);
            referenceClone = readtable(['../../Data/shahData/referenceClones/' samples{5} '_reference_clone.txt'], 'Delimiter', '\t');
            referenceClone = referenceClone.Removed_Clone + 1;
            Y = [Y(:, 1:referenceClone-1), additional_column, Y(:, referenceClone:end)];

            % Find closest indices in T for given days
            closest_indices = zeros(length(days), 1);
            for i = 1:length(days)
                [~, closest_indices(i)] = min(abs(T - days(i)));
            end

            % Index Y using closest_indices
            Y = Y(closest_indices, :);

            % Compute error
            error(p, :) = sum(sum((freqs - Y').^2, 1));
            % error(p,:) = sum(sum((freqs(:,2:end) - Y(2:end,:)').^2, 1));
        end
    else
        for p = 1:size(samples, 2)
            days = samples{2, p};
            freqs = samples{1, p};

            % Simulate data with payoff_matrix
            [T, Y] = ode45(@(t, y) replicatorEqn(t, y, payoff_matrix), days, mean(freqs(:, 1:3), 2));
            error(p, :) = sum(sum((freqs - Y').^2, 1));
            % error(p,:) = sum(sum((freqs(:,2:end) - Y(2:end,:)').^2, 1));

            % Plot the model fit
            figure();
            set(gca, 'ColorOrder', color_map, 'NextPlot', 'replacechildren');
            plot(T, Y, 'LineWidth', 4);
            hold on;

            % Plot actual data
            for i = 1:size(samples{1, p}, 1)
                plot(days, freqs(i, :), '--o', 'MarkerSize', 12, 'LineWidth', 1, ...
                    'Color', color_map(i, :), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', color_map(i, :));
            end

            % Set axis labels and title
            xlabel('Day', 'FontSize', 14, 'FontWeight', 'bold');
            ylabel('Frequency', 'FontSize', 14, 'FontWeight', 'bold');
            if ~isempty(origin_ID) && ~isempty(sample_ID)
                title([origin_ID ' ' sample_ID ' Replicate ' num2str(p) ' Model Fit'], 'FontSize', 16, 'FontWeight', 'bold');
            end

            % Add legend and grid
            legend(cellstr(samples{4, p}), 'Location', 'northeast', 'FontSize', 12);
            grid on;
            set(gca, 'GridColor', [0.9, 0.9, 0.9], 'LineWidth', 0.5);
            hold off;

            % Save the plot if OUTDIR is specified
            if ~isempty(OUTDIR)
                g = gcf;
                imageSaveName = [origin_ID '_' sample_ID '_replicate_' num2str(p) '_modelFit.png'];
                savePlace = [OUTDIR imageSaveName];
                exportgraphics(g, savePlace, 'Resolution', 300);
            end
        end
    end

    % Normalize error by the number of clones
    error = error / size(samples{1, 1}, 1);
end








