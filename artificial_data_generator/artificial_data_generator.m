%% Section 1: Generate and Save Data
clearvars
close all hidden

addpath ../replicator_eqn_optimizer/

% Number of datasets to generate
numDatasets = 20;

% Cells to store payoff matrices and ODE solutions
MCell   = cell(numDatasets,1);  % formerly pmcell
YCell   = cell(numDatasets,1);  % formerly Ycell
TCell   = cell(numDatasets,1);  % formerly Tcell

% Arrays to store meta-info about each dataset
allN     = zeros(numDatasets,1);  % formerly all_clones_gen
allOmega = zeros(numDatasets,1);  % formerly all_noise_gen
indexAtDist = zeros(numDatasets,1);  % formerly pix

distanceFromESS = cell(numDatasets,1);  % formerly d_from_ess

% Define ranges for number of clones and noise levels
minClones = 2;
maxClones = 5;
minNoise  = 0.0;
maxNoise  = 0.20;

for i = 1:numDatasets
    % Randomly pick number of clones
    n = randi([minClones, maxClones]);
    % Randomly pick noise standard deviation
    omega = (maxNoise - minNoise)*rand();

    % Generate random frequencies for initial conditions (X)
    X = rand(n,1);
    X = X / sum(X);

    % Initialize the payoff matrix (M) and validate with isESS
    validMatrix = false;  % Flag to check matrix validity

    while ~validMatrix
        M = 2*rand(n)-1;  % Generate random payoff matrix

        % Determine how many entries we want to keep
        % (At least 2n-1, but also cap at 10)
        k = min(2*n - 1, 10);

        % Flatten the matrix into a vector
        M_vec = M(:);

        % Sort the entries by absolute value in descending order
        [~, idx_sort] = sort(abs(M_vec), 'descend');

        % Identify the indices of the top k values
        top_k_idx = idx_sort(1:k);

        % Create a zero vector of the same size
        M_vec_new = zeros(size(M_vec));

        % Place the top k original values in their positions
        M_vec_new(top_k_idx) = M_vec(top_k_idx);

        % Reshape back into an n-by-n matrix
        M = reshape(M_vec_new, [n, n]);

        % If fewer than 2 non-zero entries remain, skip
        if nnz(M) < 2
            continue;
        end

        % Solve the replicator equation
        [T, Y] = ode45(@(t, y) replicatorEqn(t, y, M), [0 50], X);

        % Test for ESS
        [fit_vals, eig_vals] = isESS(M, Y(end,:)');

        % Check ESS conditions: negative real parts & fitness range <= 0.1
        if all(eig_vals(:) < 0) && range(round(fit_vals, 1)) < 0.1
            validMatrix = true;  % Matrix meets ESS conditions
        end
    end

    ess_states = findESS(M);
    d = (Y - ess_states{1,1}).^2;
    d = sqrt(mean(d,2));
    distanceFromESS{i} = d;

    [~, p_ix] = min(abs(d - 0.25));
    indexAtDist(i) = p_ix;

    % Store generated data
    MCell{i} = M;
    YCell{i} = Y;  % entire solution
    TCell{i} = T;  % timepoints

    allN(i)     = n;
    allOmega(i) = omega;
end

% Sort 'indexAtDist' in ascending order
[sortedPix, sortedIndices] = sort(indexAtDist, 'ascend');

figure()
for q = 1:numDatasets
    dee = distanceFromESS{q,1};
    imagesc('XData', q, 'YData', [0 50], 'CData', dee);
    hold on
end
hold off

% Save the generated data for later use
save('generated_data1000_11timepoints.mat', ...
     'MCell','YCell','TCell','allN','allOmega',...
     'numDatasets','distanceFromESS','-v7.3');

disp('Data generation complete. File "generated_data1000_11timepoints.mat" saved.');

%% Section 2: Load Data and Perform Inference + Violin Plots
clearvars
close all hidden

addpath ../utils/

% Load the previously generated data
load('generated_data1000_11timepoints.mat', ...
     'MCell','YCell','TCell','allN','allOmega','numDatasets','distanceFromESS');

% Pre-allocate arrays for results from the inference
MPrimeCell    = cell(numDatasets,1);   % formerly rmcell
all_timepoints= zeros(numDatasets,1);
allRMSEnonZero= zeros(numDatasets,1);  % formerly all_mean_diff
allRMSE       = zeros(numDatasets,1);  % formerly all_diff
allRMSE0      = zeros(numDatasets,1);  % formerly all_sign_diff
all_equiv     = false(numDatasets,1);

YSubCell  = cell(numDatasets,1);
TSubCell  = cell(numDatasets,1);

nonZeroCount = zeros(numDatasets,1);  % formerly nonzee
totalCount   = zeros(numDatasets,1);  % formerly numlee

% Set optimizer options
nsp = 80; % number of start points used in optimizer
opts = optimoptions(@fmincon, 'Algorithm', 'interior-point');
ms   = MultiStart('UseParallel', true, 'Display', 'off');

for i = 1:numDatasets
    M      = MCell{i};   % original payoff matrix
    Y      = YCell{i};
    T      = TCell{i};
    n      = allN(i);
    omega  = allOmega(i);

    % Identify the index up to t=30
    [~, ix] = min(abs(T - 30));

    Y = Y(1:ix,:);
    T = T(1:ix,:);

    % Generate 22 indices spanning from 1 to (up to) ix
    sampleIndices = round(linspace(1, length(T), 22));  % formerly d_ix

    % Subset T
    TSub = T(sampleIndices);
    TSubCell{i} = TSub;

    % Subset Y
    YSub = Y(sampleIndices, :);
    YSubCell{i} = YSub;

    % Use testForFreqDepEffects on the final Y (subsampled)
    [M_test, Pval, ~] = testForFreqDepEffects(TSub', YSub', [], []);
    correctedPval = Pval / size(Y,2);  % dividing by n, e.g.

    significanceScores = -log10(correctedPval) .* abs(M_test);
    [sortedScores, sortIndices] = sort(significanceScores(:), 'descend');

    % For thresholding, pick (2n-1) or simply use nnz of M
    scoreIdx = nnz(M); 
    if scoreIdx == 0
        % If M is all zeros (unlikely, but a safe check), skip
        continue;
    end
    threshold = sortedScores(scoreIdx);

    baseMatrix = (significanceScores >= threshold) .* sign(M_test);

    % Add Gaussian noise to Y (i.e., produce Y')
    YPrime = YSub + omega * randn(size(YSub));  % formerly noisy_Y
    % Ensure non-negative values
    YPrime(YPrime < 0) = 0;
    % Re-normalize each row so they sum to 1
    row_sums = sum(YPrime, 2);
    YPrime  = YPrime ./ row_sums;

    % Pack data for the inference functions
    samples = {YPrime'; TSub'; 'User-Generated-Data'; {}; 'Sample-1'; 1};

    % Bounds for the matrix entries
    ub = ones(n);
    lb = -1 * ub;

    % Find the indices we will optimize over (non-zero initial guess)
    idx_nonzero = find(baseMatrix ~= 0);

    % First optimization
    problem1 = createOptimProblem('fmincon', ...
        'objective', @(params) reo_likelihood_function(params, samples, M, []), ...
        'x0', baseMatrix(idx_nonzero), ...
        'lb', lb(idx_nonzero), 'ub', ub(idx_nonzero), ...
        'options', opts);

    rs1    = RandomStartPointSet('NumStartPoints', nsp);
    points = list(rs1, problem1);

    try
        [params1, ~] = run(ms, problem1, CustomStartPointSet(points));
        baseMatrix(idx_nonzero) = params1;

        % Re-run with updated guess: finalS
        finalS  = replicator_eqn_optimizer(samples, baseMatrix, lb, ub);
        idx2    = find(finalS~=0);
        problem2 = createOptimProblem('fmincon', ...
            'objective', @(params) reo_likelihood_function(params, samples, finalS, []), ...
            'x0', finalS(idx2), 'lb', lb(idx2), 'ub', ub(idx2), ...
            'options', opts);

        rs2     = RandomStartPointSet('NumStartPoints', nsp);
        points2 = list(rs2, problem2);
        [params2, ~] = run(ms, problem2, CustomStartPointSet(points2));

        MPrime = zeros(n);  % formerly recovered_matrix
        MPrime(idx2) = params2;

        MPrimeCell{i} = MPrime;

        % Check strategic equivalence
        [isEquivalent, resultMat] = testStrategicEquivalence(MPrime, M);
        all_equiv(i) = isEquivalent;

        % Compute errors
        new_mat = [M(:), MPrime(:)];
        % Restrict to rows where at least one is non-zero:
        new_mat = new_mat(~all(new_mat==0,2), :);

        allRMSEnonZero(i) = rmse(new_mat(:,1), new_mat(:,2));    % partial measure
        allRMSE(i)        = rmse(M(:), MPrime(:));               % full RMSE

        pm0 = (M == 0);
        rm0 = (MPrime == 0);
        allRMSE0(i) = rmse(double(pm0(:)), double(rm0(:)));  % zero-structure difference

        nonZeroCount(i) = nnz(M) + nnz(MPrime);
        totalCount(i)   = numel(M) + numel(MPrime);
        all_timepoints(i) = numel(TSub);
    catch
        warning('Dataset %d: optimizer unable to solve', i);
        continue
    end
end

disp('Optimization complete.');

%% Visualizations
% ------------------------------------------------------------------------
% Now produce the violin plots using the stored arrays
% 1) We have allRMSEnonZero, allRMSE0, all_timepoints
% 2) We also have the # of clones and noise from allN, allOmega
% ------------------------------------------------------------------------

% Define noise categories
edges = [0, 0.025, 0.1333, 0.2];
[~, noise_bins] = histc(allOmega, edges);
noise_labels    = {'Low noise','Medium noise','High noise'};
unique_clones   = unique(allN);
unique_clones   = unique_clones(unique_clones>0);

% Summarize how many non-zero entries there are
nonZeroCountsLine = [nonZeroCount, totalCount];  % formerly dynamic_error_line
keepers_idx = ~all(nonZeroCountsLine==0,2);
nonZeroCountsLine = nonZeroCountsLine(keepers_idx, :);
clone_keepers = allN(keepers_idx);

two_clone_idx   = clone_keepers==2;
p2 = mean(nonZeroCountsLine(two_clone_idx,1)./nonZeroCountsLine(two_clone_idx,2));
e2 = sqrt(2*p2*(1-p2));

three_clone_idx = clone_keepers==3;
p3 = mean(nonZeroCountsLine(three_clone_idx,1)./nonZeroCountsLine(three_clone_idx,2));
e3 = sqrt(2*p3*(1-p3));

four_clone_idx  = clone_keepers==4;
p4 = mean(nonZeroCountsLine(four_clone_idx,1)./nonZeroCountsLine(four_clone_idx,2));
e4 = sqrt(2*p4*(1-p4));

five_clone_idx  = clone_keepers==5;
p5 = mean(nonZeroCountsLine(five_clone_idx,1)./nonZeroCountsLine(five_clone_idx,2));
e5 = sqrt(2*p5*(1-p5));

% Example usage of violin plots:
createViolinPlots(allRMSEnonZero, allN, noise_bins, noise_labels, ...
    unique_clones, 'RMSE (Non-Zero Entries)', [0 1.40], sqrt(2/3));
filename = [ strrep('../../Results/RMSE_NonZero_Entries',' ','_'), '_11timepoints.png' ];
saveas(gcf, filename);

createViolinPlots(allRMSE, allN, noise_bins, noise_labels, ...
    unique_clones, 'RMSE', [0 1.40], sqrt(2/3));
filename = [ strrep('../../Results/RMSE',' ','_'), '_11timepoints.png' ];
saveas(gcf, filename);

createViolinPlots(allRMSE0, allN, noise_bins, noise_labels, ...
    unique_clones, 'Inference of Zeroes', [0 1.60], [e2,e3,e4,e5]);
filename = [ strrep('../../Results/ZeroesPlot',' ','_'), '_11timepoints.png' ];
saveas(gcf, filename);

createViolinPlots(all_timepoints, allN, noise_bins, noise_labels, ...
    unique_clones, 'Number of timepoints', [10 30], []);
filename = [ strrep('../../Results/TimepointsPlot',' ','_'), '_11timepoints.png' ];
saveas(gcf, filename);

%% ---------- Local Functions Below ---------- %%
function createViolinPlots(metric_data, all_clones, noise_bins, noise_labels, unique_clones, metric_name, ylims, cFlag)
    % Create a figure for violin plots
    figure('Color','w');
    num_clones = length(unique_clones);
    
    for i_cl = 1:num_clones
        cl = unique_clones(i_cl);
        subplot(num_clones,1,i_cl)
        hold on

        idx_cl = (all_clones == cl);

        % Extract data by noise category
        data_low = metric_data(idx_cl & noise_bins == 1);
        data_med = metric_data(idx_cl & noise_bins == 2);
        data_high= metric_data(idx_cl & noise_bins == 3);

        all_data = [data_low; data_med; data_high];
        groups_raw = [repmat({'Low noise'},   length(data_low),1);
                      repmat({'Medium noise'},length(data_med),1);
                      repmat({'High noise'},  length(data_high),1)];

        groups = categorical(groups_raw, noise_labels, 'Ordinal', true);

        % Plot the violin
        violinplot(all_data, groups, 'ShowData', true);

        % Optionally draw an expected-error line
        if ~isempty(cFlag)
            if isscalar(cFlag)
                expected_error = cFlag;
            elseif length(cFlag) == num_clones
                expected_error = cFlag(i_cl);
            else
                warning('cFlag length does not match number of unique_clones. Skipping error line.');
                expected_error = [];
            end

            if ~isempty(expected_error) && ~isnan(expected_error)
                yline(expected_error, '--', 'Expected Error', ...
                      'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
            end
        end

        title(sprintf('%d Clones', cl));
        ylabel(metric_name);
        if i_cl == num_clones
            xlabel('Noise Category');
        else
            xlabel('');
        end
        if ~isempty(ylims)
            ylim(ylims);
        end
        hold off
    end
end


function [isEquivalent, result] = testStrategicEquivalence(A, B)
    % Check if two payoff matrices A and B are strategically equivalent
    % Returns:
    %   isEquivalent (boolean)
    %   result       (chosen matrix closer to B)
    if ~isequal(size(A), size(B))
        error('Matrices A and B must have the same dimensions.');
    end

    A_flat = A(:);
    B_flat = B(:);

    % Solve B = a*A + b via least squares
    X = [A_flat, ones(size(A_flat))];
    coefficients = X \ B_flat;

    a = coefficients(1);
    b = coefficients(2);

    % Reconstruct B using the solved coefficients
    B_reconstructed = a*A + b;

    diff_B_reconstructed = sum(abs(B(:) - B_reconstructed(:)));
    diff_A = sum(abs(B(:) - A(:)));

    if diff_B_reconstructed < diff_A
        result = B_reconstructed;
        isEquivalent = false;
    else
        result = A;
        isEquivalent = true;
    end
end

