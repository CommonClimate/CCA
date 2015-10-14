# CCA
Optimized Matlab code to perform Canonical Correlation Analysis for paleoclimate reconstructions

—— cca_master_example.m: an example illustrating the workflow of climate field reconstructions using CCA. The example uses a pseudoproxy network, but the code can be generalized to use on any other datasets (including real-world proxies). The following functions are called within the master code:

—— cca_cfr.m: implement CCA reconstructions, in which the CCA parameters are estimated using the cca_cv.m function
—— cca_cfr_parallel.m: parallel the cca_cfr.m script. This is particularly useful when dealing with sparse matrices (i.e. matrices with missing values, for instance, a real-world proxy matrix). In cca_cfr_parallel.m, a set of CCA parameters is needed for each pattern of missing values. The script is coded to parallel the computation of CCA parameters of different patterns.
—— cca_cv.m: estimate the set of CCA parameters using cross-validation (CV). By default, leave-half-out CV as implemented in Smerdon et al., 2010 is used.
—— cca_bp.m: make predictions using CCA in the Barnett-Preisendorfer version (CCA-BP)
—— standardize.m: a function that standardizes a matrix and return the standardized matrix and the original matrix’s mean and standard deviation.

Reference:
Wang, J., J. Emile-Geay, D. Guillot, J.E. Smerdon, and B. Rajaratnam, 2014: Evaluating climate field reconstruction techniques using improved emulations of real-world conditions, Clim., Past, 10, 1-19.   http://www.clim-past.net/10/1/2014/cp-10-1-2014.html

Smerdon, J.E., A. Kaplan, D. Chang, and M.N. Evans, 2010:  A pseudoproxy evaluation of the CCA and RegEM methods for reconstructing climate fields of the last millennium, J. Climate, 23, 4856-4880.   http://rainbow.ldeo.columbia.edu/~alexeyk/Papers/Smerdon_etal2010CCApblshd.pdf
