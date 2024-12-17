# Dissertation
This repository contains the code files used for my dissertation.
</br>
https://nbn-resolving.org/urn:nbn:de:bsz:ch1-qucosa2-945992
</br></br>


## Chapter 3 - Direct inversion methods for the NFFT

### Visualization 
**$\bullet$ plot_nonzeros_matrixB_1d.m** </br>
$\quad\leadsto$ Code file to plot the structure of the sparse matrix $\boldsymbol B$ as done in Figure 3.1

**$\bullet$ phi.m** </br>
$\quad\leadsto$ Auxiliary code file containing the definitions of the NFFT window functions

**$\bullet$ phi_hat.m** </br>
$\quad\leadsto$ Auxiliary code file containing the definitions of the Fourier transformed window functions

**$\bullet$ plot_grids.m** </br>
$\quad\leadsto$ Code file to plot the considered grids as done in Figures 3.3 and 3.4


### Numerical examples

*Please note that the official NFFT software (https://github.com/NFFT/nfft) is necessary to run most of the following code files. </br>
A precompiled version can be downloaded from https://www-user.tu-chemnitz.de/~potts/nfft/download.php*

**$\bullet$ example_quality_densitycomp.m** </br>
$\quad\leadsto$ Code file to check the quality of the proposed density compensation method, see Figure 3.5

**$\bullet$ example_comparison_densitycomp.m** </br>
$\quad\leadsto$ Code file to compare the proposed density compensation method to methods from the literature for $d=1$, see Figures 3.6 and 3.7

**$\bullet$ example_densitycomp_phantom.m** </br>
$\quad\leadsto$ Code file to compare the proposed density compensation method to methods from the literature for $d=2$, see Figure 3.8

**$\bullet$ example_min_norm_B.m** </br>
$\quad\leadsto$ Code file to check the quality of the proposed matrix optimization method, see Figure 3.9

**$\bullet$ example_comparison_optB.m** </br>
$\quad\leadsto$ Code file to compare the proposed matrix optimization method to methods from the literature with $d=1$, see Figure 3.10

**$\bullet$ example_comparison.m** </br>
$\quad\leadsto$ Code file to compare the proposed density compensation method and the proposed matrix optimization method, see Figure 3.11



</br></br>
## Chapter 4 - Regularized Shannon sampling formulas 

### Visualization of window functions
**$\bullet$ plot_frequency_window_function.m** </br>
$\quad\leadsto$ Code file to plot the considered frequency window functions as done in Figures 4.1, 4.2 and 4.3

**$\bullet$ plot_spatial_window_function.m** </br>
$\quad\leadsto$ Code file to plot the considered spatial window functions as done in Figures 4.4, 4.5, 4.7 and 4.10

**$\bullet$ cardinal_bspline.m** </br>
$\quad\leadsto$ Auxiliary code file containing the centered cardinal B-splines 


### Visualization of properties in proofs
**$\bullet$ assertion_bspline.m** </br>
$\quad\leadsto$ Code file to visualize the assertions made in the proof for the B-spline window function, see Figure 4.6

**$\bullet$ assertion_sinh.m** </br>
$\quad\leadsto$ Code file to visualize the assertions made in the proof for the sinh-type window function, see Figures 4.8 and 4.9

**$\bullet$ assertion_cKB.m** </br>
$\quad\leadsto$ Code file to visualize the assertions made in the proof for the continuous Kaiser-Bessel window function, see Figures 4.11 and 4.12


### Numerical examples
**$\bullet$ example_nonrobustness.m** </br>
$\quad\leadsto$ Code file to demonstrate the nonbustness of the classical Shannon sampling sums, see Figure 4.13

**$\bullet$ example_psihat_lin.m** </br>
$\quad\leadsto$ Code file for the visualization of the error bound for the linear frequency window function as done in Figure 4.14

**$\bullet$ example_spatial_windows_verification.m** </br>
$\quad\leadsto$ Code file for the visualization of the error bounds for the spatial window function as done in Figures 4.16 - 4.23

**$\bullet$ example_comparison_windows.m** </br>
$\quad\leadsto$ Code file for the comparisons of frequency windows and spatial windows, see Figures 4.15, 4.24 and 4.25



</br></br>
## Chapter 5 - Fast sinc methods 

### Verification
**$\bullet$ verification_tensorization_sinc_approx.m** </br>
$\quad\leadsto$ Code file to verify the tensor product notation of the approximation of the sinc function

**$\bullet$ testfile_nfft_like_approach.m** </br>
$\quad\leadsto$ Code file to compare the NFFT-like approach and the NFFT as done in Figure 5.1


### Numerical examples
**$\bullet$ example_err_approx_sinc.m** </br>
$\quad\leadsto$ Code file to visualize the sinc approximation by means of Clenshaw-Curtis and Gauss-Legendre quadrature, see Figures 5.2 and 5.3

**$\bullet$ legpts.m** </br>
$\quad\leadsto$ Code file containing the Legendre points and weights, source: https://github.com/chebfun/chebfun

**$\bullet$ example_err_approx_sinc_least_square.m** </br>
$\quad\leadsto$ Code file to visualize the sinc approximation by means of the least squares approach, see Figure 5.4

**$\bullet$ example_maxerr_sinctrafo.m** </br>
$\quad\leadsto$ Code file to demonstrate the accuracy of the fast sinc transform, see Figure 5.5

**$\bullet$ fastsinc.m** </br>
$\quad\leadsto$ Class file containing the fast sinc transform described in Algorithm 5.10

**$\bullet$ example_sinctrafo_frequency_regularization.m** </br>
$\quad\leadsto$ Code file to demonstrate the accuracy of the fast sinc transform for frequency regularized Shannon sampling formulas, see Figure 5.6

**$\bullet$ example_sinctrafo_spatial_regularization.m** </br>
$\quad\leadsto$ Code file to demonstrate the accuracy of the fast sinc transform for spatial regularized Shannon sampling formulas, see Figure 5.7

**$\bullet$ example_approx_nfft_like.m** </br>
$\quad\leadsto$ Code file to compare the approximation error of the NFFT-like approach and the NFFT, see Figure 5.8

