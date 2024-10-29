# Random Perturbation Subsampling for Rank Regression with Massive Data

This repository contains the code used for the simulation and empirical studies presented in the paper *"Random Perturbation Subsampling for Rank Regression with Massive Data."*

## Structure

### 1. Simulation Code (`simulation code for rpsrr`)

The folder `simulation code for rpsrr` contains all the code for the simulations in the paper. It includes three subfolders corresponding to the simulation steps in the paper:

- **1 selection of m**: Code related to the selection of \(m\) in Section 4.1.
- **2 effect of stochastic weights and r**: Code for the effect of stochastic weights and subsample size \(r\) in Section 4.2.
- **3 comparison of subsampling methods**: Code for comparing different subsampling methods in Section 4.3.

### 2. Empirical Code (`real data code for dmap`)

The folder `real data code for dmap` contains the code for the empirical analysis, organized into two subfolders:

- **MSPEs**: Code for calculating the Mean Squared Prediction Errors (MSPEs), corresponding to the empirical results in Figure 5.
- **PEs**: Code for calculating the Prediction Errors (PEs), corresponding to the empirical results in Figure 6.

Additionally, the empirical data file is named `F.csv`, which supports the analysis results.

## Usage

Each folder contains code and related data files for reproducing the results presented in the paper. Please refer to the specific folder corresponding to the simulation or empirical analysis you are interested in.
