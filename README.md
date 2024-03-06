# Stage-specific-mortality-and-diversification
This repository contains the code associated with the paper "Differential stage-specific mortality as a mechanism for diversification" in The American Naturalist.

-------------------
GENERAL INFORMATION
-------------------
1. Associated article: "Differential stage-specific mortality as a mechanism for diversification", under review at The American Naturalist.


2. Author information:
   Name: Catalina Chaparro-Pedraza
   Institution: Department of Fish Ecology and Evolution, Eawag
   Email: Catalina.Chaparro@eawag.ch
   

3. Summary: MATLAB scripts used to generate the numerical results and figures presented in the article (figures 2-4).
   To produce figures 2-4 please run Make_figs.m.
   Warning: Figure 4 takes around 30 seconds for the calculations.


4. Scripts written by Catalina Chaparro-Pedraza.


5. List of files:

   - make_figs.m : This script calls functions and scripts that make the calculations (see below) and plot the results.
   - FitnessLandscapePhi.m : This script calculates the fitness gradient and curvature of the fitness function.
   - findavgmort.m : This function returns the average mortality for a given combination of demographic parameters
   - EvoDyn.m : This function returns the rate of change of the trait given the ecological parameters using eq. 4 of the paper
   - EcoEqui.m : This function encodes the ecological dynamics (eq. 1-3 in the paper) of m number of ectomorphs in an environment with n number of resources.
   - EcoEquivsPhi.m : This script calculates the population density and composition as well as the threshold of minimum productivity for diversification for a set of values of the differential mortality for diversification.
   - SimulateEcoEvoDyn.m : This script simulates the evolutionary dynamics of a population colonising an environment with 2 food resources.

More detailed information about the script can be found in the header of each file.


6. Scripts run using MATLAB R2019b (Update4, Mac)
