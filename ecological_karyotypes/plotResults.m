function [error, combinationAvgFitness] = plotResults(payoff_matrix, samples, OUTDIR, origin_ID, sample_ID, solMat)
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
        % Extract the RGB color values into color_map
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

            % Compute error between model and data
            error(p, :) = sum(sum((freqs - Y').^2, 1));
            
            % Plot the model fit (frequencies)
            figure();
            set(gca, 'ColorOrder', color_map, 'NextPlot', 'replacechildren');
            plot(T, Y, 'LineWidth', 4);
            hold on;
            for i = 1:size(samples{1, p}, 1)
                plot(days, freqs(i, :), '--o', 'MarkerSize', 12, 'LineWidth', 1, ...
                    'Color', color_map(i, :), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', color_map(i, :));
            end
            xlabel('Day', 'FontSize', 14, 'FontWeight', 'bold');
            ylabel('Frequency', 'FontSize', 14, 'FontWeight', 'bold');
            if ~isempty(origin_ID) && ~isempty(sample_ID)
                title([origin_ID ' ' sample_ID ' Replicate ' num2str(p) ' Model Fit'], ...
                      'FontSize', 16, 'FontWeight', 'bold');
            end
            legend(cellstr(samples{4, p}), 'Location', 'northeast', 'FontSize', 12);
            grid on;
            set(gca, 'GridColor', [0.9, 0.9, 0.9], 'LineWidth', 0.5);
            hold off;
            
            % Save the model fit plot if OUTDIR is specified
            if ~isempty(OUTDIR)
                g = gcf;
                imageSaveName = [origin_ID '_' sample_ID '_replicate_' num2str(p) '_modelFit.png'];
                savePlace = [OUTDIR imageSaveName];
                exportgraphics(g, savePlace, 'Resolution', 300);
            end
        end
    else
        for p = 1:size(samples, 2)
            days = samples{2, p};
            freqs = samples{1, p};

            % Simulate data with payoff_matrix
            [T, Y] = ode45(@(t, y) replicatorEqn(t, y, payoff_matrix), days, mean(freqs(:, 1:3), 2));
            error(p, :) = sum(sum((freqs - Y').^2, 1));

            % === Compute per-clone weighted fitness vectors ===
            num_clones = size(Y, 2);
            clone_fitness = zeros(size(Y, 1), num_clones);
            clone_fitness_shifted = zeros(size(Y, 1), num_clones);
            for i = 1:size(Y, 1)
                x = Y(i,:)';
                f = payoff_matrix * x;
                c = abs(min(f));
                f_shifted = f + c;
                clone_fitness(i, :) = (x .* f)';
                clone_fitness_shifted(i, :) = (x .* f_shifted)';
            end
            % Overall population fitness:
            avg_fitness = sum(clone_fitness, 2);
            avg_fitness_shifted = sum(clone_fitness_shifted, 2);
            % ===============================================

            % Plot the model fit (frequencies)
            figure();
            set(gca, 'ColorOrder', color_map, 'NextPlot', 'replacechildren');
            plot(T, Y, 'LineWidth', 4);
            hold on;
            for i = 1:size(samples{1, p}, 1)
                plot(days, freqs(i, :), '--o', 'MarkerSize', 12, 'LineWidth', 1, ...
                    'Color', color_map(i, :), 'MarkerFaceColor', 'none', 'MarkerEdgeColor', color_map(i, :));
            end
            xlabel('Day', 'FontSize', 14, 'FontWeight', 'bold');
            ylabel('Frequency', 'FontSize', 14, 'FontWeight', 'bold');
            if ~isempty(origin_ID) && ~isempty(sample_ID)
                title([origin_ID ' ' sample_ID ' Replicate ' num2str(p) ' Model Fit'], ...
                      'FontSize', 16, 'FontWeight', 'bold');
            end
            legend(cellstr(samples{4, p}), 'Location', 'northeast', 'FontSize', 12);
            grid on;
            set(gca, 'GridColor', [0.9, 0.9, 0.9], 'LineWidth', 0.5);
            hold off;

            % Save the model fit plot if OUTDIR is specified
            if ~isempty(OUTDIR)
                g = gcf;
                imageSaveName = [origin_ID '_' sample_ID '_replicate_' num2str(p) '_modelFit.png'];
                savePlace = [OUTDIR imageSaveName];
                exportgraphics(g, savePlace, 'Resolution', 300);
            end
            
            % ----- Updated: Compute population fitness for all clone combinations -----
        clones = cellstr(samples{4, p});
        % Recalculate initial condition (using the same as in the full simulation)
        init_all = mean(freqs(:, 1:3), 2);
        % Initialize containers for combination labels, average fitness values, and the fitness time series.
        combinationLabels = {};
        combinationAvgFitness = [];
        combinationFitnessTime = {};  % <-- New cell array output
        t_subs = {};
        num_clones = length(init_all);

        % Loop over combinations of clones (minimum combination size 2)
        for k = 2:num_clones
            combos = nchoosek(1:num_clones, k);
            for idx = 1:size(combos, 1)
                indices = combos(idx, :);
                % Subset the initial condition for the current clone combination.
                new_init = init_all(indices);
                new_init = new_init(:);  % Ensure a column vector for ODE45

                % Extract the corresponding submatrix of the payoff matrix.
                sub_payoff = payoff_matrix(indices, indices);

                % Re-simulate the replicator dynamics for this subset over the same days.
                [T_subset, Y_subset] = ode45(@(t, y) replicatorEqn(t, y, sub_payoff), [min(days) max(days)], new_init);

                %% Compute the fitness over time for the subset simulation
                num_timepoints_subset = size(Y_subset, 1);
                fitness_over_time = zeros(num_timepoints_subset, 1);
                shiftValue = 1e-6;  % small constant to avoid zero fitness values
                for t = 1:num_timepoints_subset
                    x_subset_norm = Y_subset(t, :);
                    fitness_over_time(t) = shiftValue + x_subset_norm * (sub_payoff * x_subset_norm');
                end

                % Store the fitness time series in the cell array.
                combinationFitnessTime{end+1} = fitness_over_time;
                t_subs{end+1} = T_subset;

                comboLabel = strjoin(clones(indices), ',');

                %% Compute average fitness over all timepoints for the bar plot
                avg_fit_subset = mean(fitness_over_time);
                combinationAvgFitness(end+1) = avg_fit_subset;
                combinationLabels{end+1} = comboLabel;
            end
        end
        end
    end

    % Normalize error by the number of clones
    error = error / size(samples{1, 1}, 1);
end

