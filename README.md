# ETKPF
The Ensemble Transform Kalman Particle Filter


## Description
Implementation of the Ensemble Transform Kalman Particle Filter, including the EnKF, ETKF and PF as special cases. The code is written as an R-package, but its main function is to quickly develop new functionalities in R and port them to fortran for used with COSMO. The documentation is rather sparse, but starting with the example in inst/examples should give an idea of the high-levels functions available. For the low-level functions, one can look at the suite of tests in tests/testthat which check that the R- and fortran versions of functions produce the same results. 

## Instruction
To install the package use devtools::install_github("robertsy/assimilr") or simply clone the repo and install the package locally. You can run the tests with devtools::test(). For every new function that you add in R and fortran you should create a test in tests/testthat/test_r_vs_f90.R

