Predictive Likelihood Clustering for Functional Data

This repository contains the R code used in the paper on predictive likelihood clustering for functional data, including simulation studies, real data applications, and plotting scripts.

Repository Contents
Simulation Examples

The repository includes representative simulation examples for the following clustering methods:

Predictive Likelihood Clustering
kmeans.fd
funFEM

To keep the repository concise, only a subset of the simulation settings reported in the manuscript are included. For example:

Predictive Likelihood: 40 curves, K = 5, sigma = 0.5
kmeans.fd: 20 curves, K = 4, sigma = 0.5
funFEM: 20 curves, K = 2, sigma = 0.5

Other simulation settings reported in the paper can be reproduced by modifying the corresponding parameters in the scripts, including:

Number of curves (N)
Number of clusters (K)
Noise level (sigma)
True cluster memberships
Data generation settings
Tuning parameters
Basis Selection Experiments

The manuscript includes experiments comparing different basis representations, including:

Fourier basis
Quadratic basis
Mixed basis
B-spline basis

The basis choice can be modified directly within the corresponding scripts to reproduce the results reported in the manuscript.

Error Distribution Experiments

The repository contains scripts corresponding to different error distributions used in the simulation studies:

Normal errors
t-distributed errors
Uniform errors

These scripts were used to generate the results reported in the manuscript.

Real Data Applications

The repository includes code for the real data examples discussed in the paper, including:

Yeast gene expression data
Vertical Density Profile (VDP) data

Plotting scripts used to generate the figures are also provided when available.

Notes

The code was developed for research purposes and reflects the workflow used during the development of the manuscript.

Some simulation results in the manuscript were obtained by repeatedly modifying parameters and rerunning the corresponding scripts. Therefore, not every individual simulation setting appearing in the manuscript is provided as a separate file. Representative examples are included, and other settings can be reproduced through parameter adjustments.

Software

The code was developed in R.

Packages used in the project may include:

cluster
fossil
splines
funFEM
fda.usc