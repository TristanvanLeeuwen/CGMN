This Matlab code can be used to reproduce the examples from the paper:

T. van Leeuwen - Fourier analysis of the CGMN method for solving the Helmholtz equation, ArXiv:1210:2644, 2012.

Contents:
* HelmND.m     - setup Helmholtz matrix
* mat2R.m      - convert Matlab sparse matrix format to band-storage format
* DKSWPR.m     - perform Kaczmarz sweeps on band-storage matrix
* sweepR_mex.c - computational kernel of DKSWPR.m
* errortest.m  - reproduces figure 1 in the paper
* cgiter       - reproduces figure 2 in the paper
* exp1         - reproduces figure 3 in the paper