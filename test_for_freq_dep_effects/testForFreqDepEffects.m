% testForFreqDepEffects.m

function [M, Pval, growthRates] = testForFreqDepEffects(days, freqs, cloneIDs, colors)
% Test for frequency-dependent effects between subpopulations in time-series data
% Input:
%   - days: Vector of time points
%   - freqs: Matrix of frequencies (rows: subpopulations, columns: time points)
%   - cloneIDs: Cell array of clone identifiers
%   - colors: Matrix of colors for plotting (rows: subpopulations, columns: RGB values)
% Output:
%   - M: Correlation matrix
%   - Pval: Matrix of p-values for the correlations

% Check inputs
if size(freqs, 2) ~= length(days)
    error('The number of columns in freqs must match the length of days.');
end
% if size(freqs, 1) ~= length(cloneIDs) || size(freqs, 1) ~= size(colors, 1)
%     error('The number of rows in freqs must match the length of cloneIDs and colors.');
% end

numSubpopulations = size(freqs, 1);
M = zeros(numSubpopulations, numSubpopulations);
Pval = zeros(numSubpopulations, numSubpopulations);

% Check if there are enough time points for meaningful analysis
if length(days) < 3
    warning('Not enough time points for meaningful analysis.');
    return;
end

% Calculate growth rates for each clone
growthRates = zeros(numSubpopulations, length(days) - 1);
for i = 1:numSubpopulations
    growthRates(i, :) = diff(log(freqs(i, :))) ./ diff(days);
end

% Check for complex numbers and replace them with zero
if any(imag(growthRates(:)) ~= 0)
    warning('Complex numbers detected in growthRates. Replacing them with zero.');
    growthRates(imag(growthRates) ~= 0) = 0;
end

% Loop over each subpopulation to investigate its fitness
for q = 1:numSubpopulations
    for v = 1:numSubpopulations
        % Calculate correlation and p-value
        [R, P] = corr(growthRates(q, :)', freqs(v, 1:end-1)', 'Type', 'Pearson');
        M(q, v) = R;  % Update the correct position in the matrix
        Pval(q, v) = P;  % Update the correct position in the matrix
    end
end

% % Plot results
% figure();
% for q = 1:numSubpopulations
%     subplot(ceil(numSubpopulations / 3), 3, q);
%     hold on;
%     for v = 1:numSubpopulations
%             scatter(freqs(v, 1:end-1), growthRates(q, :), 50, colors(v, :), 'filled', 'DisplayName', cloneIDs{v});
%     end
%     xlabel("Frequency of other Clone");
%     ylabel("Growth Rate of Clone " + cloneIDs{q});
%     legend('Location', 'best');
%     title("Growth Rate of Clone " + cloneIDs{q} + " as a Function of Other Clone Frequency");
%     hold off;
% end
end
