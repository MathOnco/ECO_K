function [Corr, Pval, growthRates] = testForFreqDepEffects(days, freqs, cloneIDs)
    % Test for frequency-dependent effects between subpopulations in time-series data
    %
    % Inputs:
    %   - days   : Vector of time points (length = number of passages)
    %   - freqs  : Matrix of frequencies (rows: subpopulations, columns: time points)
    %   - cloneIDs: Cell array of clone identifiers (not absolutely required if not used)
    %
    % Outputs:
    %   - Corr       : n-by-n correlation matrix
    %   - Pval       : n-by-n matrix of p-values
    %   - growthRates: n-by-(tau-1) matrix of growth rates

    if size(freqs, 2) ~= length(days)
        error('The number of columns in freqs must match the length of days.');
    end

    numSubpopulations = size(freqs, 1);
    Corr = zeros(numSubpopulations, numSubpopulations);
    Pval = zeros(numSubpopulations, numSubpopulations);

    if length(days) < 3
        warning('Not enough time points for meaningful analysis.');
        growthRates = [];
        return;
    end

    % ------------------------------------------------
    % Calculate growth rates for each subpopulation i
    % (alpha_{i,t} = diff(log(freq_i)) / diff(days))
    % ------------------------------------------------
    growthRates = zeros(numSubpopulations, length(days) - 1);
    for i = 1:numSubpopulations
        % log(freq(i,:)) can cause complex numbers if freq(i,:) has zeros
        % or negative values -> user must ensure data is valid
        growthRates(i, :) = diff(log(freqs(i, :))) ./ diff(days);
    end

    % Check for any imaginary parts (if freq=0 at some time, log(0) -> -Inf, etc.)
    if any(imag(growthRates(:)) ~= 0)
        warning('Complex numbers detected in growthRates. Setting them to zero.');
        growthRates(imag(growthRates) ~= 0) = 0;
    end

    % ----------------------------------
    % For each pair (q,v), compute:
    %  Corr(q,v) = correlation( alpha_{q,t}, x_{v,t} )
    % ----------------------------------
    for q = 1:numSubpopulations
        for v = 1:numSubpopulations
            [R, P] = corr( growthRates(q, :)', ...
                           freqs(v, 1:end-1)', ...
                           'Type', 'Pearson' );
            Corr(q, v) = R;
            Pval(q, v) = P;
        end
    end
end

