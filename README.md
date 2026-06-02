# README: AFT Mixture Cure Model with Competing Risks (Based on Lending Club Data)

This repository contains R code for the simulation study and real-data analysis presented in the paper. The code implements an EM algorithm for an accelerated failure time (AFT) mixture cure model with two competing events (default and prepayment), **supporting three baseline distributions**: Weibull, lognormal, or loglogistic, using interval‑censored data. The analysis focuses on Lending Club loan data.

## File structure and description

| File name | Description |
|-----------|-------------|
| `default-code.R` | **Real‑data analysis for the default event.** Preprocesses the Lending Club data, builds an interval‑censored dataset, runs an EM algorithm to estimate the mixture cure model (comparing default vs. cured), performs bootstrap inference, and plots survival curves for the default sub‑distribution and the overall survival. |
| `prepayment-code.R` | **Real‑data analysis for the prepayment event.** Similar to `default-code.R`, but removes default events and focuses on prepayment vs. cured. Includes EM estimation, bootstrap standard errors, and survival plots. |
| `example(loan).R`   | **Full real‑data example for both events simultaneously.** Constructs interval‑censored data for three event types (cured, default, prepayment), estimates the full competing‑risks mixture model using an EM algorithm, obtains bootstrap inference, and produces survival curves for each cause and the overall survival. |
| `SM-AFT (A).R`     | **Simulation study (supplementary material).** Generates synthetic data under an AFT mixture cure model with two competing risks and interval‑censoring, runs the EM algorithm on 200 simulated datasets, computes bias, MSE, coverage probabilities, and compares empirical and average standard errors. This replicates the simulation results in the paper. |