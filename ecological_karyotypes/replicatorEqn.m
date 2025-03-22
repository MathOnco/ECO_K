function xdot = replicatorEqn(~, x, payoff_matrix)
% replicatorEqn: Defines the replicator equation ODE.
% This function computes the rate of change in clone frequencies based on the replicator equation,
% as described in the Methods (see Eqn. for replicator dynamics).
    
    % Compute fitness for each clone using the payoff matrix
    % (this corresponds to f_i = sum_j x_j * M_{ij} in the Methods)
    fitness_values = payoff_matrix * x;
    % Compute the population-average fitness: \bar{f} = x' * fitness_values
    avg_fitness = x' * fitness_values;
    % Compute change in clone frequencies using the replicator equation: \dot{x}_i = x_i (f_i - \bar{f})
    xdot = x .* (fitness_values - avg_fitness);
end
