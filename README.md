
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MJMbamlss

<!-- badges: start -->
<!-- badges: end -->

The goal of MJMbamlss is to provide a working implementation of the
proposed approach found in [Volkmann, Umlauf, Greven
(2023)](https://arxiv.org/abs/2311.06409): “Flexible joint models for
multivariate longitudinal and time-to-event data using multivariate
functional principal components”.

You can find more information on the package in the README.txt file in
the inst/ folder. As a general outline for the usage of MJMbamlss refer
to the following steps:

1.  Preprocess the data to be of long format with fixed variable name
    ‘marker’ for the longitudinal outcomes factor variable.

2.  To estimate the MFPC basis, first remove observations with too
    little information. Then use wrapper function ‘preproc_MFPCA’ to
    estimate MFPCs. The number of MFPCs can be determined looking at the
    ratio of explained variance.

3.  Prepare the model formula. The formula is a list specifying each
    additive predictor separately, except for marker-specific
    predictors. That is, the baseline hazard can be specified with
    ‘Surv2(.)’ functions on the left of the ‘~’, baseline covariates are
    specified by ‘gamma ~’, error measurments with ‘sigma ~’. The alpha
    and mu predictors need to specify the model formulas with
    interactions of the variable ‘marker’, so exclude the intercept and
    specify all model terms as marker-interactions. Use the smooth terms
    ‘bs = “unc_pcre”’ for the functional principal components based
    random effects. Each ‘unc_pcre’ term needs to be supplied with an
    ‘xt’ argument ‘“mfpc”’ that contains a multiFunData object of the
    corresponding MFPC. Note also that each smooth term should contain
    the ‘xt’ argument ‘“scale” = FALSE’.

4.  Prepare the data for the model fit. Use the wrapper function
    ‘attach_wfpc’ to add evaluations of the MFPCs to the data set.

5.  Fit the model using ‘bamlss’ by specifying the family ‘mjm_bamlss’.

Please use the provided files in the folder inst/ as a reference for
your analysis.

## Installation

NOTE: For now, the package has unresolved Namespace-problems
(‘smooth.construct’ is exported as a generic function by both ‘mgcv’ and
‘bamlss’). Installation of this package via Github should work, as well
as using the package for the analysis and simulation. A fix for the
Namespace problem is under way.

<!-- You can install the stable release version of MJMbamlss from [CRAN](https://cran.r-project.org/) with: -->
<!-- ``` r -->
<!-- install.packages("MJMbamlss") -->
<!-- ``` -->

You can install the development version of MJMbamlss from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("alexvolkmann/MJMbamlss")
```
