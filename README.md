[![made-with-MATLAB](https://img.shields.io/badge/Made%20with-MATLAB_R2023b-orange)]()
[![GitHub latest commit](https://badgen.net/github/last-commit/YidanXue/stokes_2d_bifurcation)](https://GitHub.com/YidanXue/LARS/commit/)
[![arXiv](https://img.shields.io/badge/DOI-10.48550%2FarXiv.2309.11230-blue)](https://doi.org/10.48550/arXiv.2309.11230)
![X (formerly Twitter) Follow](https://img.shields.io/twitter/follow/YidanXue?link=https%3A%2F%2Ftwitter.com%2FYidanXue)

Repository for the manuscript "Stokes flows in a two-dimensional bifurcation" by Yidan Xue, Stephen J. Payne and Sarah L. Waters. The arXiv preprint can be found at https://arxiv.org/abs/2309.11230.

Prerequisites
----------------------

The codes need to be run in MATLAB or Octave, its free alternative. The recommended version is MATLAB_R2023b. The user needs to have aaa.m in the same folder (or add it in the default MATLAB path) for the AAA rational approximation. The easist way to do this is to download Chebfun at https://www.chebfun.org. 

Exmample MATLAB codes
----------------------
The codes presented in the main folder can be used to reproduce the results in bifurcations with straight channels. For the bifurcations with smooth channels, the codes can be found in the 'smooth' folder. For the bifurcations with a fixed particle, please use the codes in the 'particle' folder (except for Figure 12). The codes for the machine learning model and trained neural networks can be found in the 'ML' folder. Here is a summary of the structure of the codes and how to reproduce the figures presented in the paper.

Main folder:
1) 'bifurcation.m' computes the example simulation in Figure 2, using the lightning algorithm.
2) 'diameter_effects.m' examines the effects of channel width on flow conductance and reproduces Figure 3.
3) 'angle_cases.m' presents the flow characteristics for different bifurcation angles and reproduces Figure 4.
4) 'angle_effects.m' investigates the effects of bifurcation agnles on flow conductance and reproduces Figure 5.
5) 'pressure_cases.m' visualises the flow characteristics for different outlet pressures (Figure 11).
6) 'pressure_cases_particle.m' investigates the same scenarios but with a fixed cylindrical particle in the centre (Figure 12).

The 'smooth' folder:
1) 'smooth.example.m' reproduces the example simulation in Figure 6, using the AAA algorithm.
2) 'diameter_effects_smooth.m' reproduces the effects of channel width on flow conductance of a bifurcation with curved boundaries (Figure 7).

The 'particle' folder:
1) 'particle.m' computes the example case in Figure 8, using the lightning algorithm with a series method.
2) 'cylinder_location.m' investigates the particle location on flow conductance (Figure 9).

The 'ML' folder:
There are 3 subfolders, where 'bifurcation', 'smooth_boundary' and 'particle' reproduce Figures 10a, 10b and 10c, respectively. In each folder, 'parameters.mat' stores the input parameters for both training and validation. Running 'data_generation.m' will simulate the flow conductance for each input parameter set using the LARS algorithm, and output 'g0c.mat', 'g1c.mat' and 'g2c.mat'. Finally, executing 'ml_fitting.m' will train the 3 neural networks 'g0_ml.mat', 'g1_ml.mat' and 'g2_ml.mat' using the first 80\% of the simulation results, and compute the loss using the last 20\% of the simulation results.
