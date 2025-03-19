% shahVignette

clearvars
close all hidden
cd '~/Repositories/Frequency-dependentSelection_paper/replicator_eqn_optimizer/'
import bioma.data.*
addpath ../test_for_freq_dep_effects/
addpath ../bootstrapping/
addpath ../
OUTDIR = 'Results/';

% Initialize cell arrays for storing payoff matrices and ESS results
pmcell = {};
ecell = {};

fprintf('Script started running at: %s\n', datetime('now'));

% Start a parallel pool for optimization (as many optimization steps run in parallel)
if isempty(gcp('nocreate'))
    parpool('local'); % Starts a parallel pool with the default number of workers
end

% Get main directory containing datasets; each "origin" corresponds to a replicate group
dirinfo = dir('paths/*origin*');

% Initialize table to summarize fit results across datasets
fitSummaryTable = table('Size', [0, 20], 'VariableTypes', ...
    {'cell', 'cell', 'cell', 'cell', 'double', 'double', 'double', ...
    'double', 'double', ...
    'double', 'double', 'double', 'double', ...
    'double', 'double', 'double', 'double', 'double','double','double'}, 'VariableNames', ...
    {'datasetName', 'label', 'sampleName', 'repID', 'error', 'AIC', 'BIC', ...
    'matEntryAbsMean', 'cloneCount', ...
    'interactionsDetected', 'timepoints', 'rxFracs', 'matEntryParams', ...
    'matEntryZeros', 'matEntryMax', 'matEntryMin', 'matEntryMean', ...
    'ns','ns2','ns3'});

% Loop through each replicate group (origin)
for v = 1:length(dirinfo)
    % Get local directory containing samples for the current origin
    localPath = ['paths/' dirinfo(v).name '/*sample*'];
    localDir = dir(localPath);

    % Loop through each sample within the current origin
    for l = 1:length(localDir)
        datasetName = {};
        label = {};
        sampleName = {};
        repID = {};
        rxFracs = {};
        samples = {};
        solMat = {};

        % Define path to the replicate directory
        replicatePath = ['paths/' dirinfo(v).name '/' localDir(l).name];
        if isfolder(replicatePath)
            replicateDir = dir(replicatePath);
            replicateDir = replicateDir(3:end,:); % Skip . and ..
            if strcmp(replicateDir(1).name, '.DS_Store')
                replicateDir(1) = [];
            end

            % Identify excluded clones based on frequency threshold (10%) across all replicates,
            % as described in the methods (a clone never exceeding 10% frequency is excluded).
            for m = 1:length(replicateDir)
                p = [replicatePath '/' replicateDir(m).name];
                inputTable = readtable(p, 'ReadRowNames', 0);
                inputTable = removevars(inputTable, {'Var1', 'replicateID', 'rxFrac'});
                if m == 1
                    exclude = find(all(table2array(inputTable) < 0.10, 2));
                else
                    exclude = intersect(exclude, find(all(table2array(inputTable) < 0.10, 2)));
                end
            end

            % Process each replicate within the current sample
            for m = 1:length(replicateDir)
                p = [replicatePath '/' replicateDir(m).name];
                inputTable = readtable(p, 'ReadRowNames', 0);
                cloneIDs = table2array(inputTable(:, 1));
                replicateID = inputTable.replicateID(1);
                rxFracs{m} = inputTable.rxFrac(1);
                inputTable = removevars(inputTable, {'Var1', 'replicateID', 'rxFrac'});
                % Exclude clones that never exceed 10% frequency (see methods section)
                keepers = setdiff(1:size(inputTable, 1), exclude);
                inputTable = inputTable(keepers, :);
                cloneIDs = cloneIDs(keepers);
                sampleID = extractBefore(replicateDir(m).name, '_');
                % Define time vector (days) based on passages (30-day intervals)
                days = 30 * (1:size(inputTable, 2)) - 30;
                % Normalize frequencies so that each column sums to 1
                X = table2array(inputTable);
                X = X ./ sum(X);
                % Store labels and dataset identifiers for later use
                label{m} = strjoin(unique(cellfun(@(x) extractAfter(x, "_"), inputTable.Properties.VariableNames, 'UniformOutput', false), 'stable'), '');
                datasetName{m} = [dirinfo(v).name '_' localDir(l).name '_replicate_' num2str(m)];
                repID{m} = [dirinfo(v).name '_' localDir(l).name];
                sampleName{m} = sampleID;
                samples = [samples {X; days; sampleID; cloneIDs; repID{m}; l}];
            end

            % Initialize variables to store metrics from the optimization and model fitting
            matEntryAbsMean = nan(size(samples,2),1);
            timepoints = nan(size(samples,2),1);
            matEntryParams = nan(size(samples,2),1);
            interactionsDetected = nan(size(samples,2),1);
            matEntryZeros = nan(size(samples,2),1);
            matEntryMax = nan(size(samples,2),1);
            matEntryMin = nan(size(samples,2),1);
            matEntryMean = nan(size(samples,2),1);
            cloneCount = nan(size(samples,2),1);

            % Test for frequency-dependent effects between pairs of clones.
            % This step computes significance scores for pairwise interactions (see Methods: "Identifying which pairs of clones are candidates for frequency-dependent interactions")
            n = size(samples{1, 1}, 1);
            if length(days)>3
                ems = zeros(n, n);
                pvals = zeros(n, n);
                for p = 1:size(samples, 2)
                    [Corr, Pval, growthRates] = testForFreqDepEffects(samples{2, p}, samples{1, p}, samples{4, p});
                    pvals = pvals + Pval;
                    ems = ems + Corr;
                end
                Pval = pvals / size(samples, 2);
                % Pval = Pval / size(samples, 2);
                Corr = ems / size(samples, 2);

                numTests        = n * n; 
                correctedPval    = Pval * numTests;
                % correctedPval   = min(rawBonfPvals, 1);  

                % Calculate significance scores based on correlation and p-value (see methods section)
                significanceScores = -log10(correctedPval) .* abs(Corr);
                [sortedScores, sortIndices] = sort(significanceScores(:), 'descend');
                [sortedQ, sortedV] = ind2sub(size(significanceScores), sortIndices);
                sortedInteractions = table(sortedQ, sortedV, sortedScores, ...
                    'VariableNames', {'Subpopulation1', 'Subpopulation2', 'SignificanceScore'});

                score_idx = (2*n)-1;

                % Initialize payoff matrix M_initial based on significance scores.
                % This represents the initial non-zero interactions to be optimized (see Methods for the payoff matrix M)
                M_initial = significanceScores >= max(sortedScores(score_idx));
                M_initial = M_initial .* sign(Corr);

                % Visualize the significance of interactions using a heatmap
                figure()
                heatmap(significanceScores);
                colormap('summer');
                sgtitle('Significance of Interactions');
                g = gcf;
                imageSaveName = [dirinfo(v).name '_' localDir(l).name '_interactionsHeatmap.png'];
                savePlace = [OUTDIR imageSaveName];
                exportgraphics(g, savePlace, 'Resolution', 300);
            else
                % If insufficient time points, default to identity matrix for M_initial
                M_initial=zeros(n);
                M_initial(1:size(M_initial,1)+1:end) = 1;
            end

            % Set bounds for the payoff matrix parameters, reflecting that fitness interactions are bounded between -1 and 1.
            ub = ones(n);
            lb =-1*ub;

            % Create initial optimization problem using fmincon with the SQP algorithm.
            % This step corresponds to estimating the payoff parameters in the replicator equation (Methods: optimization routine)
            nsp = 80; % number of start points used in optimizer problem
            opts = optimoptions(@fmincon, 'Algorithm', 'sqp');
            ms = MultiStart('UseParallel', true);
            ms.Display = 'off';

            idx = find(M_initial~=0);

            problem = createOptimProblem('fmincon', 'objective', ...
                @(params) reo_likelihood_function(params, samples, M_initial, []), ...
                'x0', M_initial(idx), 'lb', lb(idx), 'ub', ub(idx), 'options', opts);
            rs = RandomStartPointSet('NumStartPoints', nsp);
            points = list(rs, problem);

            % Run initial optimization to estimate nonzero payoff parameters by minimizing negative log-likelihood
            [params, negative_log_likelihood] = run(ms, problem, CustomStartPointSet(points));

            M_initial(idx) = params;

            % Further refine the payoff matrix using the replicator_eqn_optimizer function,
            % which iteratively removes non-significant interactions based on AIC/BIC (see Methods)
            finalM = replicator_eqn_optimizer(samples, M_initial, lb, ub); %, dirinfo(v).name, localDir(l).name);

            idx = find(finalM~=0);

            % Update bounds based on standard deviation of optimized parameters
            sd=std(finalM(idx));
            ub = finalM+sd';
            lb = finalM-sd';

            % Re-run optimization with updated bounds to further refine parameter estimates
            problem = createOptimProblem('fmincon', 'objective', ...
                @(params) reo_likelihood_function(params, samples, finalM, []), ...
                'x0', finalM(idx), 'lb', lb(idx), 'ub', ub(idx), 'options', opts);
            rs = RandomStartPointSet('NumStartPoints', nsp);
            points = list(rs, problem);

            [params, negative_log_likelihood] = run(ms, problem, CustomStartPointSet(points));

            % Calculate model selection criteria (AICc and BIC) using the optimized likelihood,
            % as described in the methods section.
            [aicc, bic] = reo_calculate_AICBIC(negative_log_likelihood, params, [samples{2,:}]);

            % Construct the final payoff matrix with the optimized parameters.
            payoff_matrix = zeros(n);
            payoff_matrix(idx) = params;

            % Store the final payoff matrix and evolutionary stable strategy (ESS) for each sample
            pmcell{v,l} = payoff_matrix;
            ecell{v,l}={findESS(payoff_matrix)};

            % Plot the dynamics of the payoff matrix, corresponding to the replicator dynamics simulation (Methods)
            figure()
            matrixDynamicsPlot(payoff_matrix, samples{4, 1},dirinfo(v).name)
            title([dirinfo(v).name localDir(l).name "Matrix Dynamics"]);
            fontsize(16, "points");
            g = gcf;
            imageSaveName = [dirinfo(v).name '_' localDir(l).name '_matrixDynamicsPlot.png'];
            savePlace = [OUTDIR imageSaveName];
            exportgraphics(g, savePlace, 'Resolution', 300);

            % Plot solution of the replicator equation and compute error between predicted and observed frequencies
            error = reo_plotResults(payoff_matrix, samples, OUTDIR, dirinfo(v).name, localDir(l).name, []);

            % Run bootstrapping to estimate confidence in the interactions (model robustness),
            % as described in the methods section.
            num_bootstrap = 50;
            [bootstrap_estimates, num_significant,num_significant2,num_significant3] = reo_bootstrap(payoff_matrix, samples, num_bootstrap, ub, lb);

            % Gather and store additional output results for each replicate
            AIC = nan(size(samples,2),1);
            BIC = nan(size(samples,2),1);
            ns = nan(size(samples,2),1);
            ns2 = nan(size(samples,2),1);
            ns3 = nan(size(samples,2),1);
            for p = 1:size(samples,2)
                matEntryAbsMean(p,:) = mean(abs(finalM(finalM~=0)));
                timepoints(p,:) = size(samples{1, p}, 2);
                matEntryParams(p,:) = length(finalM(finalM~=0));
                interactionsDetected(p,:) = length((finalM(finalM~=0)));
                matEntryZeros(p,:) = length((finalM(finalM==0)));
                matEntryMax(p,:) = max(finalM(:));
                matEntryMin(p,:) = min(finalM(:));
                matEntryMean(p,:) = mean(finalM(finalM~=0));
                cloneCount(p,:) = size(samples{1, p}, 1);
                AIC(p,:) = aicc;
                BIC(p,:) = bic;
                ns(p,:) = num_significant;
                ns2(p,:) = num_significant2;
                ns3(p,:) = num_significant3;
            end

            % Compile results into a table for export
            output = table(error, AIC, BIC, matEntryAbsMean, cloneCount, interactionsDetected, ...
                timepoints, matEntryParams, matEntryZeros, matEntryMax, ...
                matEntryMin, matEntryMean, ns, ns2, ns3);

            datasetName = datasetName';
            label = label';
            sampleName = sampleName';
            repID = repID';
            rxFracs=rxFracs';
            fst_lead = table(datasetName, label, sampleName, repID, rxFracs);
            fst = [fst_lead output];
            fitSummaryTable = [fitSummaryTable; fst];
        end
    end
end

% Write the compiled summary results to an Excel file
writetable(fitSummaryTable, 'shahDataFitSummary.xlsx')

% Close the parallel pool
delete(gcp('nocreate'));

