
# Contents -------------------------------------------------------
# The two folders contain the code to reproduce the simulation and the analysis
# of the PBC data set.


# Application -------------------------------
# pbc_application.R: Prepare the PBC data, estimate MFPCs, fit the MJM with
#   bamlss and JMbayes2.


# Simulation --------------------------------
# Two folders 'scenarioI' and 'scenarioII' corresponding to the described
# simulation scenarios with files:
# *_data.R: Data generation of simulation
# *_bamlss_tru.R: Estimate the MJM with bamlss using the data generating MFPCs
# *_bamlss_est1.R: Estimate the MJM with bamlss using estimated MFPCs without
#                  truncation
# *_bamlss_est99.R: Estimate the MJM with bamlss using estimated MFPCs with
#                  truncation
# *_jmbayes.R: Estimate the MJM with JMbayes2
# sim_eval.R: Evaluate the simulation results




# Model fitting outline ------------------------------------------
# In order to fit the proposed MJM approach, follow these steps:

# 1. ----------------------------------------
# Preprocess the data to be of long format with fixed variable name 'marker' for
# the longitudinal outcomes factor variable).

# 2. ----------------------------------------
# To estimate the MFPC basis, first remove observations with too little
# information. Then use wrapper function 'preproc_MFPCA' to estimate MFPCs. The
# number of MFPCs can be determined looking at the ratio of explained variance.

# 3. ----------------------------------------
# Prepare the model formula. The formula is a list specifying each additive
# predictor separately, except for marker-specific predictors. That is, 
# the baseline hazard can be specified with 'Surv2(.)' functions on the left of
# the '~', baseline covariates are specified by 'gamma ~', error measurments
# with 'sigma ~'. The alpha and mu predictors need to specify the model 
# formulas with interactions of the variable 'marker', so exclude the intercept
# and specify all model terms as marker-interactions. Use the smooth terms
# 'bs = "unc_pcre"' for the functional principal components based random 
# effects. Each 'unc_pcre' term needs to be supplied with an 'xt' argument
# '"mfpc"' that contains a multiFunData object of the corresponding MFPC. Note
# also that each smooth term should contain the 'xt' argument '"scale" = FALSE'.

# 4. ----------------------------------------
# Prepare the data for the model fit. Use the wrapper function 'attach_wfpc' to
# add evaluations of the MFPCs to the data set.

# 5. ----------------------------------------
# Fit the model using 'bamlss' by specifying the family 'mjm_bamlss'.
