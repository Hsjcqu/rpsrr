# README

This repository contains the R code used for the simulation and application in the paper *"Random Perturbation Subsampling for Rank Regression with Massive Data"*. The code is organized into two main directories: one for simulation studies and one for application.

## Repository Structure

### 1. R Code for Simulation

This folder includes the R scripts used in the simulation studies discussed in the thesis. The folder is structured as follows:

*   **function8.R**:\
    Contains the primary functions used in the simulations described in the thesis.

*   **simulationforMm.R**:\
    Code for the simulation presented in Section 4.1, "*The Effect of The Selection of m*."

*   **simulationforv.R**:\
    Code for the simulation presented in Section 4.2, "*The Effect of Stochastic Weights and Subsample Size r*."

*   **Comparison/**:\
    Contains scripts for simulations in Section 4.3, "*Comparison with Several Existing Subsampling Methods*," including:

    *   **simulation11-13.R**: Simulates Algorithm 2 under Cases 1-3.

    *   **simulationforA.R**: Simulates Algorithm 4 under Cases 1-3.

    *   **simulationfortime.R**: Calculates CPU time for the proposed algorithm.

### 2. R Code for Real Data

This folder includes the R scripts used in the empirical analysis section of the thesis.

## How to Use

### Prerequisites

*   Ensure that R is installed on your machine.

*   Additional R packages may be required based on the specific scripts. Refer to the comments within each script for any package dependencies.

### Running Simulations

1.  Navigate to the `R code for simulation` folder.

2.  To replicate results from Section 4.1, execute the `simulationforMm.R` script.

3.  Similarly, run other scripts for their corresponding sections.

### Application (Real Data Analysis)

1.  Navigate to the `R code for real data` folder.

2.  Run the provided scripts to replicate the empirical analysis.

