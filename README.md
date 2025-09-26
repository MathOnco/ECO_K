# ECO_K

ECO_K is a comprehensive MATLAB toolkit designed to infer frequency-dependent clonal interaction networks from longitudinal population data. It employs an evolutionary game theory framework, modeling clonal dynamics using the replicator equation. The primary output is a **payoff matrix** that quantifies how different karyotype-defined subpopulations promote or inhibit each other's growth, revealing the underlying competitive and cooperative relationships.

The toolkit is particularly suited for analyzing time-series data from biological systems like cancer cell populations, microbial communities, or other evolving systems where clonal frequency is tracked over time.



---

## Core Method & Workflow

The main script, `shahVignette.m`, is a vignette that implements the analytical pipeline described in:

> Salehi, S., Kabeer, F., Ceglia, N. et al. Clonal fitness inferred from time-series modelling of single-cell cancer genomes. Nature 595, 585–590 (2021). https://doi.org/10.1038/s41586-021-03648-3

The automated workflow proceeds through three main stages:

1.  **Initial Interaction Screen**: A correlation-based test (`testForFreqDepEffects`) identifies a set of plausible interactions to initialize the model.

2.  **Beam Search Model Selection**: An efficient search algorithm (`ecological_karyotypes`) systematically simplifies the model to find the most parsimonious interaction network with the best Bayesian Information Criterion (BIC).

3.  **Bootstrapping**: Statistical bootstrapping (`bootstrap_func`) is used to calculate confidence intervals and assess the significance of the final inferred interaction parameters.



---

## Getting Started

### 1. Dependencies

* **MATLAB** (R2022b or later)
* **MATLAB Toolboxes**:
    * Parallel Computing Toolbox
    * Optimization Toolbox
    * Statistics and Machine Learning Toolbox

### 2. Data Structure

Place your frequency data in a nested directory structure within the `paths/` folder: `paths/{origin}/{sample}/{replicate}.csv`.

### 3. Execution

With your data in place, run the main script from the MATLAB command window:

```matlab
shahVignette
