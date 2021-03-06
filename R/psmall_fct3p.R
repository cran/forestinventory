

psmall_fct3p<- function(formula.s0, formula.s1, data, phase_id, cluster, small_area, boundary_weights, exhaustive, progressbar, psmall){

  # errorchecking for psmall:
  if(small_area[["unbiased"]]==FALSE){stop("'unbiased' cannot be set to 'FALSE' for 'psmall' = 'TRUE'")}
  if(is.na(small_area[["sa.col"]])){stop("'psmall' can only be applied for objects of class 'smallarea'")}

  ## Change unbiased into FALSE for applying snthetic estimator:
  small_area[["unbiased"]] <- FALSE

  ## Call sysnthetic estimator:
  unextended <- small_area_looper_3p(formula.s0, formula.s1, data, phase_id, cluster, small_area, boundary_weights, exhaustive, progressbar, psmall)

  ## Change Point estimate of snthetic estimation into small area point estimation:
  unextended[["estimation"]][["estimate"]] <- unextended[["estimation"]][["estimate"]] + unextended[["mean_Rc_x_hat_G"]][["mean_full_Rc_x_hat_G"]]


  ## Extract everything from the output that we need to calculate the variance of the small area estimator:

  # --> for term 1:
  Z_bar_1G.list <- as.list(as.data.frame(t(unextended[["Z_bar_1G"]][,-1])))
  coefficients.list <- as.list(as.data.frame(t(unextended[["coefficients"]][["alpha"]][,-1])))
  n2<- unique(unextended[["estimation"]][["n2"]])
  n1<- unique(unextended[["estimation"]][["n1"]])

  # --> for term 2:
  Z_bar_0G.list <- as.list(as.data.frame(t(unextended[["Z_bar_0G"]][,-1])))

  # --> for term 4:
  R_bar.list <- as.list(unextended[["mean_Rc_x_hat_G"]][["mean_full_Rc_x_hat_G"]])
  Mx_s2G_bar.list <- lapply(unextended[["Mx_s2G"]], FUN=mean)
  n2G.list <- as.list(unextended[["estimation"]][["n2G"]])


  ## warning messages:
  if (any(unlist(n2G.list) == 1)) {warning("n2G for at least one Small Area is 1. Numerical calculation of 'psmall' not possible, set to NaN")}
  if (any(unlist(n2G.list) == 0)) {warning("n2G for at least one Small Area is 0. 'psmall' - estimation not possible, set to 'NA'")}



  ## calculate variance of small area estimator:

  # term1:
  term1<- mapply(function(Z, SB){t(as.matrix(Z)) %*% as.matrix(SB[["cov_beta_s2"]]) %*% as.matrix(Z)}, Z=Z_bar_1G.list, SB=unextended[["cov_coef"]])
  term1<- (1 - (n2/n1)) * term1

  ## term2:
  term2<-  mapply(function(Z, SB){t(as.matrix(Z)) %*% as.matrix(SB[["cov_alpha_s2"]]) %*% as.matrix(Z)}, Z=Z_bar_0G.list, SB=unextended[["cov_coef"]])
  term2<- (n2/n1) * term2

  # term3: only exists in the non-exhaustive case (i.e. if 'exhaustive'=NA):
  if(all(is.na(exhaustive))){
    term3 <- mapply(function(B, SZ){t(as.matrix(B)) %*% as.matrix(SZ) %*% as.matrix(B)}, B=coefficients.list, SZ=unextended[["cov_Z1_bar_0G"]])
  } else {
    term3<- 0
  }

  # term4:
  term4 <- mapply(function(n, R, R_bar, M, M_bar){    (1/(n*(n-1))) * sum( (M/M_bar)^2 * (R[["resid_full_G"]] - R_bar)^2 ) },
                  n=n2G.list,
                  R=unextended[["Rc_x_hat_G"]],
                  R_bar=R_bar.list,
                  M=unextended[["Mx_s2G"]],
                  M_bar=Mx_s2G_bar.list)


  ## Change variance of synthetic estimation into variance of small area estimator:
  unextended[["estimation"]][["g_variance"]] <- term1 + term2 + term3 + term4
  result<- unextended

  ## Change Method tp "psmall":
  result[["input"]][["method"]]<- "psmall"

  return(result)

}
