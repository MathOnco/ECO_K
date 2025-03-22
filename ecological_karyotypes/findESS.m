function [NE_list, support_list] = findESS(A, tol)
%FINDALLSYMMETRICNE  Find all symmetric mixed Nash equilibria for an n x n payoff matrix.
% 
%   [NE_list, support_list] = findAllSymmetricNE(A, tol)
%
%   INPUTS:
%       A   : (n x n) payoff matrix for a symmetric game with n pure strategies.
%             Rows and columns both index the same set of strategies.
%       tol : (optional) numerical tolerance (defaults to 1e-12).
%
%   OUTPUTS:
%       NE_list     : cell array of row‐vectors, each an n‐dim mixed strategy.
%       support_list: cell array of corresponding supports (index sets).
%
%   METHOD:
%   1) Enumerate all nonempty subsets S of {1,...,n}.
%   2) Impose "equal payoffs" among the strategies in S:
%         payoff(i) = payoff(j) for all i,j in S
%      plus the probability simplex condition:
%         sum_{i in S} x_i = 1,  x_i >= 0,  x_i=0 for i not in S.
%   3) Solve for x_S, discard any negative solutions, and check that no out‐of‐support
%      strategy has strictly greater payoff than those in S.
%   4) Collect all distinct solutions (within a tolerance).

    if nargin < 2
        tol = 1e-12;
    end

    n = size(A,1);
    if size(A,2) ~= n
        error('Matrix A must be square (n x n).');
    end

    % We will return all solutions in cell arrays
    NE_list = {};
    support_list = {};

    % Loop over all 2^n - 1 nonempty subsets of {1,...,n}
    for subsetMask = 1:(2^n - 1)
        % Extract which strategies are in the support S
        S = find(bitget(subsetMask,1:n));  % e.g. if subsetMask=5 => binary '0101' => S={1,3}

        % Number of strategies in this support
        m = length(S);

        %------------------------------------------------------------------
        % Build the linear system that forces payoffs of all in-support
        % strategies to be equal to payoff(S(1)).
        %
        % Let x(S) = unknown probabilities for the strategies in S.
        % The payoff to playing pure strategy i is P(i) = sum_{j in S} A(i,j)*x_j.
        % We want:
        %   P(S(k)) - P(S(1)) = 0   for k = 2..m
        % plus the constraint sum_{j in S} x_j = 1.
        %------------------------------------------------------------------
        M = zeros(m, m);
        b = zeros(m, 1);

        for eq = 1:(m-1)
            iRow = S(eq+1);
            iRef = S(1);
            % Each column col in 1..m represents x_{S(col)} in the system
            for col = 1:m
                jStrat = S(col);
                M(eq,col) = A(iRow,jStrat) - A(iRef,jStrat);
            end
        end
        % Probability‐sum row: sum_{j in S} x_j = 1
        M(m,:) = 1; 
        b(m) = 1;

        % Solve M*xS = b (least‐squares or direct if nonsingular).
        % Then we will check feasibility below.
        xS_candidate = M \ b;  % might be over/under‐determined, but we check validity next.

        % Negative or large rounding error => discard
        if any(xS_candidate < -tol)
            continue;
        end
        % Force small negative rounding to 0
        xS_candidate(xS_candidate<0 & xS_candidate>-tol) = 0;

        % Renormalize in case sum is slightly off
        ssum = sum(xS_candidate);
        if abs(ssum - 1) > 1e-6
            continue;  % not valid
        end
        xS_candidate = xS_candidate / ssum;

        % Build the full n‐vector x, zero outside S
        x_candidate = zeros(1,n);
        for idx = 1:m
            x_candidate(S(idx)) = xS_candidate(idx);
        end

        %------------------------------------------------------------------
        % Check that no out‐of‐support strategy has strictly higher payoff
        % than those used in S.
        % Let commonPayoff = payoff to S(1).
        % For k in outOfS, must have payoff(k) <= commonPayoff.
        %------------------------------------------------------------------
        outS = setdiff(1:n, S);
        commonPayoff = A(S(1),:) * x_candidate.';  % payoff(S(1)) against x
        isValid = true;
        for k = outS
            payoff_k = A(k,:) * x_candidate.';
            if payoff_k > commonPayoff + tol
                isValid = false;
                break;
            end
        end

        if ~isValid
            continue;
        end

        %------------------------------------------------------------------
        % We have found a candidate equilibrium.  Check if we already stored
        % it (within tolerance).  If not, add it.
        %------------------------------------------------------------------
        alreadyListed = false;
        for ii=1:length(NE_list)
            if norm(NE_list{ii} - x_candidate,1) < 1e-8
                alreadyListed = true;
                break;
            end
        end

        if ~alreadyListed
            NE_list{end+1} = x_candidate; %#ok<AGROW>
            support_list{end+1} = S;       %#ok<AGROW>
        end
    end
end
