# This function helps the "three phase" estimator function when the first phase sample size is finite.
# This corresponds to the paper on "two-phase sampling with partially exhaustive information" in global estimation
# The external and g-weight variances are calulated

# INPUTS:
#   formula.s0 = formula for "reduced" model calculable for all of s0
#   formula.s1 = formula for "full" model calculable for all of s1
#   data = data.frame
#   boundary_weights = vector of weights representing the proportion of the forest area in the s1 interpretation area (vector of 1s if NA)
#
# OUTPUTS:
#   estimate = point estimate
#   ext_variance = external variance
#   g_variance = g-weight variance
#   n0 = large ausxiliary sample size (--> reduced model)
#   n1 = small auxiliary sample size  (--> full model)
#   n2 = terrestrial sample size
#   r.squared_reduced = R squared of reduced model defined by formula.s0
#   r.squared_full = R squared of full model defined by formula.s1


global_nonexhaustive3p <- function(formula.s0, formula.s1, data, phase_id, boundary_weights,...){

  # -------------------------------------------------------------------------- #

  # retrieve phase.columnname and indicator of s2 grid id and terrestrial-grid id:
  phase.col<- phase_id[["phase.col"]]
  s1.ind<- phase_id[["s1.id"]]
  s2.ind<- phase_id[["terrgrid.id"]]

  # -------------------------------------------------------------------------- #

  # fit models with lm()-function, based on terrestrial sample s2:
  model.object.reduced<- lm(formula.s0, data=data[data[,phase.col] == s2.ind,], x=TRUE, y=TRUE) # "reduced" model
  model.object.full<-    lm(formula.s1, data=data[data[,phase.col] == s2.ind,], x=TRUE, y=TRUE) # "full" model

  # -------------------------------------------------------------------------- #

  # derive design-matrices restricted to terrestrial sample s2:
  design_matrix_1.s2<- model.object.reduced$x  #"reduced" design matrix defined at s2
  design_matrix.s2<-   model.object.full$x     #"full" design matrix defined at s2

  # derive design-matrices restricted to reduced aux.sample s1:
  design_matrix_1.s1<- design_matrix.s1_return(formula=formula.s0, data=data[data[,phase.col] %in% c(s1.ind, s2.ind),]) #"reduced" design matrix defined at s1
  design_matrix.s1<-   design_matrix.s1_return(formula=formula.s1, data=data[data[,phase.col] %in% c(s1.ind, s2.ind),]) #"full" design matrix defined at s1

  # derive design-matrices restricted to reduced aux.sample s0 (only reduced model available):
  design_matrix_1.s0 <- design_matrix.s1_return(formula=formula.s0, data=data) #"reduced" design matrix defined at s0

  # -------------------------------------------------------------------------- #

  # determine sample sizes:
  n0 <- nrow(design_matrix_1.s0)
  n1 <- nrow(design_matrix.s1)
  n2 <- nrow(design_matrix.s2)

  # -------------------------------------------------------------------------- #

  # calculate the Z_bars...:
  Z_bar_1_s0 <- apply(design_matrix_1.s0, 2, mean)
  Z1_bar_s1<- apply(design_matrix_1.s1, 2, mean) # new
  Z_bar_s1<-  apply(design_matrix.s1, 2, mean)   # new


  # calculate the Z_bars as weighted auxiliary means if boundary_weights are used:
  if(!is.na(boundary_weights)){

    Z_bar_1_s0<- boundaryweight_fct_2p3p(formula=formula.s0,
                                         model.object=model.object.reduced,
                                         data_select_from=data,
                                         boundary_weights)

    Z1_bar_s1<- boundaryweight_fct_2p3p(formula=formula.s0,
                                        model.object=model.object.reduced,
                                        data_select_from=data[data[,phase.col] %in% c(s1.ind, s2.ind),],
                                        boundary_weights)

    Z_bar_s1<- boundaryweight_fct_2p3p(formula=formula.s1,
                                       model.object=model.object.full,
                                       data_select_from=data[data[,phase.col] %in% c(s1.ind, s2.ind),],
                                       boundary_weights)
  }


  # calculate covariance matrix of Z_bar_1_s0: (page 10, [27])
  Z_1_s0_centered <- t(apply(design_matrix_1.s0, 1, function(Zx,...){Zx - Z_bar_1_s0}))
  cov_Z_bar_1_s0 <- (t(Z_1_s0_centered) %*% Z_1_s0_centered)/(n0*(n0-1))# eq [64] p. 26 Mandallaz 2p/3p technical report

  # -------------------------------------------------------------------------- #

  # calculate the g-weights:
  A1_s1_inv<- solve( (t(design_matrix_1.s1)%*%design_matrix_1.s1) / n1 )
  A_s2_inv<-  solve( (t(design_matrix.s2)%*%design_matrix.s2) / n2 )
  g<-    Z_bar_s1 %*% A_s2_inv %*% t(design_matrix.s2)
  g_1<-  Z_bar_1_s0 %*% A1_s1_inv %*% t(design_matrix_1.s2) # only needed for s2-sample


  # get regression coefficients:
  alpha<- model.object.reduced$coefficients
  beta<-  model.object.full$coefficients

  # get residuals:
  resid_reduced<- model.object.reduced$residuals # that's R_1 in the report
  resid_full<-    model.object.full$residuals    # that's R in the report


  # asymp. consistent design-based covariance-matrices  of regression coefficients:
  cov_beta_s2<-   A_s2_inv %*% ((t(resid_full*design_matrix.s2) %*% (resid_full*design_matrix.s2)) / n2^2) %*% A_s2_inv
  cov_alpha_s2<-  A1_s1_inv %*% ((t(resid_reduced*design_matrix_1.s2) %*%  (resid_reduced*design_matrix_1.s2)) / n2^2) %*% A1_s1_inv


  # -------------------------------------------------------------------------- #

  # calculate estimates:
  estimate<-     ((Z_bar_1_s0 - Z1_bar_s1) %*% alpha) + (Z_bar_s1 %*% beta)

  ext_variance<- (1/n0)*var(design_matrix_1.s0 %*% alpha) + (1/n1)*var(resid_reduced) +
                 (1-(n2/n1))*(1/n2)*var(resid_full)

  g_variance<- (alpha %*% cov_Z_bar_1_s0 %*% alpha) + (1/(n1*n2)) * sum(g_1^2 *resid_reduced^2) +
               (1 - (n2/n1)) * (1/n2^2) *sum(g^2 * resid_full^2)

#   # alternative version of design-based variance using the covariance matrices:
#   g_variance <- (alpha %*% cov_Z_bar_1_s0 %*% alpha) + ((n2/n1)*(t(Z_bar_1_s0) %*% cov_alpha_s2 %*% Z_bar_1_s0)) +
#                     ((1-(n2/n1))*(t(Z_bar_s1) %*% cov_beta_s2 %*% Z_bar_s1))



  ## ------- create outputs ------------------------------------------------- ##

  # summarize sample size info:
  samplesizes<- data.frame(cbind (n0, n1, n2))
  colnames(samplesizes)<- c("n0", "n1", "n2")
  rownames(samplesizes)<- "plots"

  estimation<- data.frame(estimate=estimate, ext_variance=ext_variance, g_variance=g_variance,
                 n0=samplesizes$n0, n1=samplesizes$n1, n2=samplesizes$n2,
                 r.squared_reduced=summary(model.object.reduced)$r.squared,
                 r.squared_full=summary(model.object.full)$r.squared)

  # ... to store inputs used:
  inputs<- list()
  inputs[["data"]]<- data
  inputs[["formula.s0"]]<- formula.s0
  inputs[["formula.s1"]]<- formula.s1
  inputs[["boundary_weights"]]<- boundary_weights
  inputs[["method"]]<- "non-exhaustive"
  inputs[["cluster"]]<- FALSE
  inputs[["exhaustive"]]<- FALSE

  # save warning-messages:
  warn.messages<- NA

  result<- list(input=inputs,
                estimation=estimation,
                samplesizes=samplesizes,
                coefficients=list(alpha=alpha, beta=beta),
                cov_alpha_s2=cov_alpha_s2,
                cov_beta_s2=cov_beta_s2,
                Z_bar_1_s0=Z_bar_1_s0,
                Z1_bar_s1=Z1_bar_s1,
                Z_bar_s1=Z_bar_s1,
                cov_Z_bar_1_s0=cov_Z_bar_1_s0,
                resid_reduced=resid_reduced,
                resid_full=resid_full,
                warn.messages=warn.messages)

  class(result)<-c("global", "threephase")

  result

}



