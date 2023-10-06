# Estimation of Stochastic Differential Equations mixed models using the Delta method approximation

R code to estimate model parameters of the Stochastic Gompertz model, considering:
 1) The fixed effects closed formula
 2) The delta approximation method for random alpha case
 3) The delta approximation method for random beta case 
 4) The delta approximation method both alpha and beta random

These codes are the support codes for the papers cited in References. 
- Model 1) follows equation (8) of reference [1].
- Model 2) follows equation (24) of reference [3].
- Model 3) follows equation (38) of reference [3].
- Model 4) follows equation (18) of reference [1].

# Instructions for running the code

1) Run the R function DA_method_function.R
2) The function has the inputs Da.fit(Mat_Time, Mat_Obs,  random, c(a,t,b,o,s))
where:
- Mat_Time - matrix of m observation time instants. Each column represents the observation time instants of an individual.
- Mat_Obs - matrix of m individual growth measure. Each column represents  the observations of an individual.
  Notes: Mat_Time and Mat_Obs must have the same dimensions. If different animals had different observation times,  the remain observations of a given column should have NA.

- random  - the random effects in the drift. If random=0, no random effects. If random=1,  a random effect on alpha
                                           If random=11,  a random effect on beta. If random=2,  a random effect on both alpha and beta

Notes about the starting parameter values for the minimization:

  - c(a,t,b,o,s): initial values for a - alpha, t - theta (give the value 0 if no random effect is used)
                                     b - beta, o - omega (give the value 0 if no random effect is used)
                                     s - sigma (the diffusion coefficient assumed fixed
 - If random=0 use c(a,b,s). if random=1 use c(a,t,b,s). if random=11 use c(a,b,o,s). if random=2 use c(a,t,b,o,s).  

# Example
We propose a data file for example. The file have simulated information on the weigth of 30 animals, taken at different age instants.


# References
 - [1] Jamba NT, Jacinto G, Filipe PA, Braumann CA. Likelihood Function through the Delta Approximation in Mixed SDE Models. Mathematics. 2022; 10(3):385. https://doi.org/10.3390/math10030385
 - [2] Jamba NT, Jacinto G, Filipe PA, Braumann CA. Stochastic differential equations mixed model for individual growth with the inclusion of genetic characteristics. Submited.
 - [3] Jamba NT, Filipe PA, Jacinto G,  Braumann CA. Estimation for stochastic differential equations mixed models using approximation methods. Submited.
