% shahVignette

clearvars
close all hidden
cd '~/Repositories/ECO_K/ecological_karyotypes/'
import bioma.data.*
addpath ../test_for_freq_dep_effects/
addpath ../bootstrapping/
addpath ../
OUTDIR = 'Results/';

pmcell = {};
ecell = {};

fprintf('Script started running at: %s\n', datetime('now'));

% Start a parallel pool (adjust the number of workers as needed)
if isempty(gcp('nocreate'))
    parpool('local'); % Starts a parallel pool with the default number of workers
end

% Get main directory
dirinfo = dir('paths/*origin*');

fitSummaryTable = table('Size', [0, 22], 'VariableTypes', ...
    {'cell', 'cell', 'cell', 'cell', 'double', 'double', 'double', ...
    'double', 'double', ...
    'double', 'double', 'double', 'double', ...
    'double', 'double', 'double', 'double', 'double', 'double', 'double','double','double'}, 'VariableNames', ...
    {'datasetName', 'label', 'sampleName', 'repID', 'error', 'AIC', 'BIC', ...
    'matEntryAbsMean', 'cloneCount', ...
    'interactionsDetected', 'timepoints', 'rxFracs', 'matEntryParams', ...
    'matEntryZeros', 'matEntryMax', 'matEntryMin', 'matEntryMean', ...
    'fitclone_aicc', 'fitclone_bic', 'ns','ns2','ns3'});


% Loop through each origin
for v = 1:length(dirinfo)
    % Get local directory
    localPath = ['paths/' dirinfo(v).name '/*sample*'];
    localDir = dir(localPath);

    % Loop through each sample in the origin
    for l = 1:length(localDir)
        datasetName = {};
        label = {};
        sampleName = {};
        repID = {};
        rxFracs = {};
        samples = {};
        solMat = {};

        % Get replicate directory
        replicatePath = ['paths/' dirinfo(v).name '/' localDir(l).name];
        if isfolder(replicatePath)
            replicateDir = dir(replicatePath);
            replicateDir = replicateDir(3:end,:); % Skip . and ..
            if strcmp(replicateDir(1).name, '.DS_Store')
                replicateDir(1) = [];
            end

            % Identify excluded clones
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

            % Process replicates
            for m = 1:length(replicateDir)
                p = [replicatePath '/' replicateDir(m).name];
                if v ==1 || v==3
                    modified_path = strrep(p, 'paths', 'fitCloneOutput');
                    modified_path = strrep(modified_path, 'cloneFreqs.csv', 'infer_x.tsv');
                    data_table = readtable(modified_path, 'FileType', 'text', 'Delimiter', '\t');
                    data_table = removevars(data_table, {'Var1', 'np'});
                    % Multiply the entries in the 'time' column by 500
                    data_table.time = data_table.time * 500;
                    % Extract unique time points
                    time_vector = unique(data_table.time);
                    % Extract unique clone IDs
                    cloneIDs = unique(data_table.K);
                    % Initialize the solution matrix
                    solution_matrix = zeros(length(time_vector), length(cloneIDs));
                    % Fill the solution matrix with clonal frequencies
                    for i = 1:length(cloneIDs)
                        cloneID = cloneIDs(i);
                        for j = 1:length(time_vector)
                            time_point = time_vector(j);
                            % Find the clonal frequency for the current clone at the current time point
                            idx = find(data_table.K == cloneID & data_table.time == time_point);
                            if ~isempty(idx)
                                % Assuming there is only one matching row, pick the first one
                                solution_matrix(j, i) = data_table.X(idx(1));
                            end
                        end
                    end
                    solMat=[solMat {solution_matrix; time_vector}];
                end
                inputTable = readtable(p, 'ReadRowNames', 0);
                cloneIDs = table2array(inputTable(:, 1));
                replicateID = inputTable.replicateID(1);
                rxFracs{m} = inputTable.rxFrac(1);
                inputTable = removevars(inputTable, {'Var1', 'replicateID', 'rxFrac'});
                keepers = setdiff(1:size(inputTable, 1), exclude);
                inputTable = inputTable(keepers, :);
                cloneIDs = cloneIDs(keepers);
                sampleID = extractBefore(replicateDir(m).name, '_');
                days = 30 * (1:size(inputTable, 2)) - 30;
                freqs = table2array(inputTable);
                freqs = freqs ./ sum(freqs);
                label{m} = strjoin(unique(cellfun(@(x) extractAfter(x, "_"), inputTable.Properties.VariableNames, 'UniformOutput', false), 'stable'), '');
                datasetName{m} = [dirinfo(v).name '_' localDir(l).name '_replicate_' num2str(m)];
                repID{m} = [dirinfo(v).name '_' localDir(l).name];
                sampleName{m} = sampleID;
                samples = [samples {freqs; days; sampleID; cloneIDs; repID{m}; l}];
            end

            fitclone_aicc = NaN(size(samples,2),1);
            fitclone_bic = NaN(size(samples,2),1);
            if ~isempty(solMat)
                fcLikelihood=likelihood_function(0,samples,0,solMat);
                for p=1:size(samples,2)
                    [fitclone_aicc(p), fitclone_bic(p)] = calculate_AICBIC(fcLikelihood, ones(2*size(samples{1, 1}, 1)-1,1), samples{2,p});
                end
            end

            matEntryAbsMean = nan(size(samples,2),1);
            timepoints = nan(size(samples,2),1);
            matEntryParams = nan(size(samples,2),1);
            interactionsDetected = nan(size(samples,2),1);
            matEntryZeros = nan(size(samples,2),1);
            matEntryMax = nan(size(samples,2),1);
            matEntryMin = nan(size(samples,2),1);
            matEntryMean = nan(size(samples,2),1);
            cloneCount = nan(size(samples,2),1);

            % Test for frequency-dependent effects between pairs of clones
            n = size(samples{1, 1}, 1);
            colors = linspecer(n);
            if length(days)>3
                numReplicates = size(samples, 2);  % Number of replicates

                % Initialize with the first replicate
                combinedDays = samples{2,1};
                combinedFreqs = samples{1,1};

                if numReplicates>1

                    for nr = 2:numReplicates
                        currentDays = samples{2,nr};
                        currentFreqs = samples{1,nr};

                        % Find new timepoints in the current replicate that are not in combinedDays
                        newIdx = ~ismember(currentDays, combinedDays);

                        % Append these new timepoints and the corresponding frequency columns
                        combinedDays = [combinedDays, currentDays(newIdx)];
                        combinedFreqs = [combinedFreqs, currentFreqs(:, newIdx)];
                    end

                end

                % Now run the test on the combined data
                % (Assuming cloneIDs remain the same across replicates; you can use samples{4,1} for cloneIDs)
                [M, Pval, growthRates] = testForFreqDepEffects(combinedDays, combinedFreqs, samples{4,1}, colors);

                correctedPval = Pval;
                significanceScores = -log10(correctedPval) .* abs(M);
                [sortedScores, sortIndices] = sort(significanceScores(:), 'descend');
                [sortedQ, sortedV] = ind2sub(size(significanceScores), sortIndices);
                sortedInteractions = table(sortedQ, sortedV, sortedScores, ...
                    'VariableNames', {'Subpopulation1', 'Subpopulation2', 'SignificanceScore'});

                score_idx = (2*n)-1;

                basematrix = significanceScores >= max(sortedScores(score_idx));
                basematrix = basematrix .* sign(M);

                figure()
                heatmap(significanceScores);
                colormap('summer');
                sgtitle('Significance of Interactions');
                g = gcf;
                imageSaveName = [dirinfo(v).name '_' localDir(l).name '_interactionsHeatmap.png'];
                savePlace = [OUTDIR imageSaveName];
                exportgraphics(g, savePlace, 'Resolution', 300);
            else
                basematrix=zeros(n);
                basematrix(1:size(basematrix,1)+1:end) = 1;
            end

            lb_M = -1 * ones(n); % Lower bound for M (interaction matrix)
            ub_M = 1 * ones(n); % Upper bound for M

            idx_m = find(basematrix ~= 0);

            x0 = basematrix(idx_m); % Use M values and best r from previous step
            lb = lb_M(idx_m);
            ub = ub_M(idx_m);

            % Create optimization problem and set up all options
            nsp = 80; % number of start points used in optimizer problem
            opts = optimoptions(@fmincon, 'Algorithm', 'sqp');
            ms = MultiStart('UseParallel', true);
            ms.Display = 'off';

            idx = find(basematrix~=0);

            problem = createOptimProblem('fmincon', 'objective', ...
                @(params) likelihood_function(params, samples, basematrix, []), ...
                'x0', x0, 'lb', lb, 'ub', ub, 'options', opts);
            
            rs = RandomStartPointSet('NumStartPoints', nsp);
            points = list(rs, problem);

            [params, negative_log_likelihood] = run(ms, problem, CustomStartPointSet(points));

            basematrix(idx) = params;

            % Example Call
            beam_width = 3; 
            [finalM, finalBIC] = ecological_karyotypes(samples, basematrix, lb_M, ub_M, beam_width);

            finalS = finalM;
            idx = find(finalS~=0);

            sd=std(finalS(idx));
            ub_M = finalS+sd';
            lb_M = finalS-sd';

            x0 = finalS(idx); % Use M values from previous step
            lb = lb_M(idx);
            ub = ub_M(idx);

            problem = createOptimProblem('fmincon', 'objective', ...
                @(params) likelihood_function(params, samples, finalS, []), ...
                'x0', x0, 'lb', lb, 'ub', ub, 'options', opts);
            nsp=400;
            rs = RandomStartPointSet('NumStartPoints', nsp);
            points = list(rs, problem);

            [params, negative_log_likelihood] = run(ms, problem, CustomStartPointSet(points));

            [aicc, bic] = calculate_AICBIC(negative_log_likelihood, params, [samples{2,:}]);


            % check out final matrix
            payoff_matrix = zeros(n);
            payoff_matrix(idx) = params;

            figure()
            matrixDynamicsPlot(payoff_matrix, samples{4, 1},dirinfo(v).name)
            title([dirinfo(v).name localDir(l).name "Matrix Dynamics"]);
            fontsize(16, "points");
            g = gcf;
            imageSaveName = [dirinfo(v).name '_' localDir(l).name '_matrixDynamicsPlot.png'];
            savePlace = [OUTDIR imageSaveName];
            exportgraphics(g, savePlace, 'Resolution', 300);

            pmcell{v,l} = payoff_matrix;

            % plot solution and get error for each replicate
            [error, combinationAvgFitness] = plotResults(payoff_matrix, samples, OUTDIR, dirinfo(v).name, localDir(l).name, []);

            % Run bootstrapping
            num_bootstrap = 1000;

            ub = ub_M;
            lb = lb_M;

            [bootstrap_estimates, num_significant,num_significant2,num_significant3] = bootstrap_func(payoff_matrix, samples, num_bootstrap, ub, lb);

            % gather up other output results
            AIC = nan(size(samples,2),1);
            BIC = nan(size(samples,2),1);
            ns = nan(size(samples,2),1);
            ns2 = nan(size(samples,2),1);
            ns3 = nan(size(samples,2),1);
            for p = 1:size(samples,2)
                matEntryAbsMean(p,:) = mean(abs(finalS(finalS~=0)));
                timepoints(p,:) = size(samples{1, p}, 2);
                matEntryParams(p,:) = length(finalS(finalS~=0));
                interactionsDetected(p,:) = length((finalS(finalS~=0)));
                matEntryZeros(p,:) = length((finalS(finalS==0)));
                matEntryMax(p,:) = max(finalS(:));
                matEntryMin(p,:) = min(finalS(:));
                matEntryMean(p,:) = mean(finalS(finalS~=0));
                cloneCount(p,:) = size(samples{1, p}, 1);
                AIC(p,:) = aicc;
                BIC(p,:) = bic;
                ns(p,:) = num_significant;
                ns2(p,:) = num_significant2;
                ns3(p,:) = num_significant3;
            end

            output = table(error, AIC, BIC, matEntryAbsMean, cloneCount, interactionsDetected, ...
                timepoints, matEntryParams, matEntryZeros, matEntryMax, ...
                matEntryMin, matEntryMean, fitclone_aicc, fitclone_bic, ns, ns2, ns3);

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

    writetable(fitSummaryTable, 'Results/fitSummaryTable.xlsx')

    delete(gcp('nocreate'));

end
