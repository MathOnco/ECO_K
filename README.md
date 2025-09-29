# ECO_K

```ECO_K``` is a comprehensive MATLAB toolkit designed to infer frequency-dependent clonal interaction networks from longitudinal population data. It employs an evolutionary game theory framework, modeling clonal dynamics using the replicator equation. The primary output is a **payoff matrix** that quantifies how different karyotype-defined subpopulations promote or inhibit each other's growth, revealing the underlying competitive and cooperative relationships.

The toolkit is particularly suited for analyzing time-series data from biological systems like cancer cell populations, microbial communities, or other evolving systems where clonal frequency is tracked over time.

---

## Example: Fitting the Hawk-Dove Game

This example shows how to use the core functions of the repository to infer a known payoff matrix from synthetically generated data. It serves as a minimal working example to illustrate the underlying logic. You can run this code in a single script once the repository's functions are on your MATLAB path.

### 1. Generate Synthetic Data from a Known Model

First, we'll define the classic Hawk-Dove payoff matrix, generate an ideal frequency trajectory using the replicator equation, and add some noise to simulate a real dataset.

```matlab
% --- Define the Hawk-Dove Game and Generate Data ---
% V=Value of resource (2), C=Cost of fighting (10). Assume C > V.
payoff_matrix_true = [ (2-10)/2 , 2 ; 
                       0       , 2/2 ];  % True Matrix A = [-4, 2; 0, 1]

% Solve for the frequency trajectory over 20 time steps
t_span = 0:20;
[t, freqs] = ode45(@(t,x) replicatorEqn(t, x, payoff_matrix_true), t_span, [0.6; 0.4]);

% Transpose and add noise to simulate experimental data
freqs = freqs'; % Clones in rows, time in columns
freqs_noisy = freqs + 0.01 * randn(size(freqs));
freqs_noisy = max(min(freqs_noisy, 1), 0);         % Clamp between 0 and 1
freqs_noisy = freqs_noisy ./ sum(freqs_noisy, 1);  % Re-normalize

% Package data into the required 'samples' cell array format
samples = {freqs_noisy; t_span'; 'Hawk-Dove'; {'Hawk', 'Dove'}}; % {frequencies; days; origin; clone IDs}
```
---

### 2. Infer the Payoff Matrix from Data

Next, we'll use the noisy data to infer the four parameters of the 2x2 payoff matrix. We do this by defining an objective function that calls likelihood_function and minimizing it with fmincon.

```matlab
% The ecological_karyotypes function requires an initial matrix, parameter bounds,
% and a beam width for the search.

% Define lower and upper bounds for the payoff matrix elements.
lower_bounds = -1 * ones(2);
upper_bounds =  ones(2);

% To get a good starting point (M_initial), we first generate a data-driven
% guess using testForFreqDepEffects.
[M_guess, ~, ~] = testForFreqDepEffects(t_span, freqs_noisy, {'Hawk', 'Dove'}, []);
initial_guess = M_guess(:); % Vectorize the matrix to use as a starting point

% We then refine this guess with a quick optimization to get a high-quality
% starting matrix for the beam search.
structure_matrix = ones(2);
objective_fun = @(params) likelihood_function(params, samples, structure_matrix, []);
preliminary_params = fmincon(objective_fun, initial_guess, [], [], [], [], lower_bounds(:), upper_bounds(:));
M_initial = reshape(preliminary_params, 2, 2);

% Set the beam width for the search algorithm.
beam_width = 3;

% Run the main function. This may take a moment as it is a comprehensive search.
% It will start with the M_initial model and test if removing any of the four
% parameters results in a model with a better (lower) BIC score.
[payoff_matrix_inferred, best_bic] = ecological_karyotypes(samples, M_initial, lower_bounds, upper_bounds, beam_width);

fprintf('Model selection complete. Best BIC found: %f\n', best_bic);
```
---

### 3. Results
Finally, plot the fit to data and display the inferred matrix alongside the true matrix to verify that the fitting procedure successfully recovered the game's parameters.

```matlab

plotResults(payoff_matrix_inferred, samples, [], [], [], []);

% --- Display the True and Inferred Matrices ---
disp('True Payoff Matrix:');
disp(payoff_matrix_true);

disp('Inferred Payoff Matrix:');
disp(payoff_matrix_inferred);
```
---

### Core Functions
This repository contains a full pipeline for automated analysis, but the Hawk-Dove example above relies on these core components:

```testForFreqDepEffects.m```: Performs a preliminary statistical screen by correlating each clone's growth rate with every other clone's frequency to generate a data-driven initial guess for the interaction matrix.

```replicatorEqn.m```: Defines the ordinary differential equation (ODE) for the replicator dynamics. This function takes the current population frequencies and a payoff matrix and returns the rate of change for each clone's frequency, forming the core of the evolutionary model.

```likelihood_function.m```: The objective function for optimization. It uses an ODE solver (which calls ```replicatorEqn.m```) to predict frequency dynamics and computes the model's negative log-likelihood against the observed data.

```ecological_karyotypes.m```: Implements a Beam Search algorithm to perform model selection and find the interaction matrix with the best BIC score for more complex models.

```bootstrap_func.m```: Manages bootstrapping to assess the statistical confidence of inferred interaction parameters.

---

### For Complex Datasets
The example above is ideal for understanding the core fitting mechanism. For analyzing complex, multi-sample experimental datasets, the main automated script ```shahVignette.m``` is recommended. It handles data loading, preprocessing, model selection, and bootstrapping across many samples as described in: 

> Salehi, S., Kabeer, F., Ceglia, N. et al. Clonal fitness inferred from time-series modelling of single-cell cancer genomes. Nature 595, 585–590 (2021). https://doi.org/10.1038/s41586-021-03648-3

---

### Dependencies
MATLAB (R2022b or later)

MATLAB Toolboxes:

Parallel Computing Toolbox

Optimization Toolbox

Statistics and Machine Learning Toolbox

