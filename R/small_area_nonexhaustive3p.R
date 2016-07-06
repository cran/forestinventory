



small_area_nonexhaustive3p<- function(formula.s0, formula.s1, data, phase_id, small_area, boundary_weights, thresh.n2g){



  # intialize warning-message vector (i.e. collector of warnings produced during function execution):
  w1<- NA; w2<- NA; w3<- NA

  # -------------------------------------------------------------------------- #

  # retrieve phase.columnname and indicator of s1 grid id and terrestrial-grid id and smallarea-columname:
  phase.col<- phase_id[["phase.col"]]
  sa.col<- small_area[["sa.col"]]
  s1.ind<- phase_id[["s1.id"]]        # s1.id identifies the small sample of the auxiliary vars (s1-sample, refers to formula.s1, ie. the full-model)
  s2.ind<- phase_id[["terrgrid.id"]]  # identifies the terrestrial sample (s2-sample)

  # get number of plots for phases in F (entire domain):
  n2<- sum(data[[phase.col]] == s2.ind)
  n1<- sum(data[,phase.col] %in% c(s1.ind, s2.ind))
  n0<- nrow(data)

  # -------------------------------------------------------------------------- #

  # fit models with lm()-function, based on entire terrestrial sample s2:
  model.object.reduced<- lm(formula.s0, data=data[data[,phase.col] == s2.ind,], x=TRUE, y=TRUE)     # "reduced" model
  model.object.full<-    lm(formula.s1, data=data[data[,phase.col] == s2.ind,], x=TRUE, y=TRUE)     # "full" model

  # extract r-squares for both models:
  r.squared_reduced<- summary(model.object.reduced)$r.squared
  r.squared_full<-    summary(model.object.full)$r.squared

  # get global residuals:
  resid_reduced<- model.object.reduced$residuals # that's R_1 in the report
  resid_full<-    model.object.full$residuals    # that's R in the report

  # get regression coefficients:
  alpha<- model.object.reduced$coefficients
  beta<-  model.object.full$coefficients

  # save indices of s2G-points in model.object.xxx:
  s2G.in.modobject.ind<- which(data[data[,phase.col] == s2.ind, ] [[sa.col]] == 1)


  # -------------------------------------------------------------------------- #
  # derive design-matrices:

  # derive design-matrices restricted to terrestrial sample s2 for entire domain F:
  design_matrix_1.s2<- model.object.reduced$x  #"reduced" design matrix defined at terrestrial points s2
  design_matrix.s2<-   model.object.full$x     #"full" design matrix defined at terrestrial points s2


  # derive design-matrices restricted to reduced aux.sample s1 or entire domain F:
  design_matrix_1.s1<- design_matrix.s1_return(formula=formula.s0, data=data[data[,phase.col] %in% c(s1.ind, s2.ind),]) #"reduced" design matrix defined at s1


  # derive design-matrices restricted to terrestrial sample s1 for small area G:
  design_matrix_1.s1G<- design_matrix.s1_return(formula=formula.s0, data=data[ data[,phase.col] %in% c(s1.ind, s2.ind) & data[,sa.col] == 1, ])
  design_matrix.s1G<-   design_matrix.s1_return(formula=formula.s1, data=data[ data[,phase.col] %in% c(s1.ind, s2.ind) & data[,sa.col] == 1, ])


  # derive design-matrices restricted to reduced aux.sample s0 (only reduced model available):
  design_matrix_1.s0G <- design_matrix.s1_return(formula=formula.s0, data=data[data[,sa.col] == 1, ]) #"reduced" design matrix defined at s0


  # determine sample sizes:
  n2G<-  sum(data[,phase.col] %in% s2.ind & data[,sa.col] == 1)
  n1G<-  nrow(design_matrix_1.s1G)
  n0G<-  nrow(design_matrix_1.s0G)

  # -------------------------------------------------------------------------- #

  # calculate the Z_bars...:
  Z_bar_1_s0G <- apply(design_matrix_1.s0G, 2, mean)
  Z1_bar_s1G<-   apply(design_matrix_1.s1G, 2, mean)
  Z_bar_s1G<-    apply(design_matrix.s1G, 2, mean)


  # calculate the Z_bars as weighted auxiliary means if boundary_weights are used:
  if(!is.na(boundary_weights)){

    Z_bar_1_s0G<- boundaryweight_fct_2p3p(formula=formula.s0,
                                         model.object=model.object.reduced,
                                         data_select_from=data[data[,sa.col] == 1, ],
                                         boundary_weights)

    Z1_bar_s1G<- boundaryweight_fct_2p3p(formula=formula.s0,
                                        model.object=model.object.reduced,
                                        data_select_from=data[ data[,phase.col] %in% c(s1.ind, s2.ind) & data[,sa.col] == 1, ],
                                        boundary_weights)

    Z_bar_s1G<- boundaryweight_fct_2p3p(formula=formula.s1,
                                        model.object=model.object.full,
                                        data_select_from=data[ data[,phase.col] %in% c(s1.ind, s2.ind) & data[,sa.col] == 1, ],
                                        boundary_weights)
  }



  # calculate covariance matrix of Z_bar_1_s0: (page 10, [27])
  Z_1_s0G_centered <- t(apply(design_matrix_1.s0G, 1, function(Zx,...){Zx - Z_bar_1_s0G}))
  cov_Z_bar_1_s0G <- (t(Z_1_s0G_centered) %*% Z_1_s0G_centered)/(n0G*(n0G-1))# eq [64] p. 26 Mandallaz 2p/3p technical report

  # -------------------------------------------------------------------------- #

  # calculate the A's: --> depends on cases "no. n2G in G"


  # ------------------------------------ #
  # 1) case "no n2g in G":

  if(n2G==0){

    if(!small_area[["unbiased"]]) {
      # if "unbiased" == FALSE (i.e. if synthetic est. applied):

      # calculate the A's: --> we have no indicator introduced and hence we can invert the design-matrix although n2G == 0:
      A1_s1_inv<- solve( (t(design_matrix_1.s1)%*%design_matrix_1.s1) / n1 )
      A_s2_inv<-  solve( (t(design_matrix.s2)%*%design_matrix.s2) / n2 )

      # get residuals in small area:
      resid_reduced_G<- NA
      resid_full_G<- NA

    } else { # if "unbiased" == TRUE:

      # calculate the A's: not possible, since we cannot invert the design-matrix:
      A1_s1_inv<- matrix(NA, nrow=ncol(design_matrix_1.s1), ncol=ncol(design_matrix_1.s1))
      A_s2_inv<-  matrix(NA, nrow=ncol(design_matrix.s2), ncol=ncol(design_matrix.s2))

      # set residuals in small area to NA:
      resid_reduced_G<- NA
      resid_full_G<- NA

      # set r-squares to NA to force entire output to appear as 'NA':
      r.squared_reduced<- NA
      r.squared_full<- NA

      w1<- warning(paste("Unbiased estimation for Small Area", small_area[["areas"]],"not possible, set to 'NA': Small Area",small_area[["areas"]],"does not contain any terrestrial data"), call. = F)

    }
  }

  # ------------------------------------ #
  # 2) case "only 1 n2g in G":

  if(n2G==1){

    # calculate the A's:
    A1_s1_inv<- solve( (t(design_matrix_1.s1)%*%design_matrix_1.s1) / n1 )
    A_s2_inv<-  solve( (t(design_matrix.s2)%*%design_matrix.s2) / n2 )

    # get residuals in small area:
    resid_reduced_G<- model.object.reduced$residuals[s2G.in.modobject.ind]
    resid_full_G<- model.object.full$residuals[s2G.in.modobject.ind]

    w2<- warning(paste("Terrestrial sample size in small area",small_area[["areas"]],"is 1. Computation of external variance not possible"), call. = F)

  }

  # ------------------------------------ #
  # 3) "normal" case (few n2g available):

  if(n2G>=2){

    # calculate the A's:
    A1_s1_inv<- solve( (t(design_matrix_1.s1)%*%design_matrix_1.s1) / n1 )
    A_s2_inv<-  solve( (t(design_matrix.s2)%*%design_matrix.s2) / n2 )

    # get residuals in small area:
    resid_reduced_G<- model.object.reduced$residuals[s2G.in.modobject.ind]
    resid_full_G<- model.object.full$residuals[s2G.in.modobject.ind]

  }

  # ------------------------------------ #


  # -------------------------------------------------------------------------- #

  # asymp. consistent design-based covariance-matrices  of regression coefficients:
  cov_beta_s2<-   A_s2_inv %*% ((t(resid_full*design_matrix.s2) %*% (resid_full*design_matrix.s2)) / n2^2) %*% A_s2_inv
  cov_alpha_s2<-  A1_s1_inv %*% ((t(resid_reduced*design_matrix_1.s2) %*%  (resid_reduced*design_matrix_1.s2)) / n2^2) %*% A1_s1_inv

  # -------------------------------------------------------------------------- #

  # calculate estimates:
  estimate<-  ((Z_bar_1_s0G - Z1_bar_s1G) %*% alpha) + (Z_bar_s1G %*% beta)

  g_variance<- (alpha %*% cov_Z_bar_1_s0G %*% alpha) +
               (n2/n1) * (t(Z_bar_1_s0G) %*% cov_alpha_s2 %*% Z_bar_1_s0G) +
               (1 - (n2/n1)) * (t(Z_bar_s1G) %*% cov_beta_s2 %*% Z_bar_s1G)


  # external variance:
  ext_variance<- (1 / n0G) * var( design_matrix_1.s0G  %*% alpha ) +
                 (1 / n1G) * var(resid_reduced_G) +
                 (1 - (n2G/n1G)) * (1/n2G) * var(resid_full_G)


  ## ------- create outputs ------------------------------------------------- ##


  # summarize sample size info:
  samplesizes<- data.frame(cbind (n0G, n1G, n2G, n0, n1, n2))
  colnames(samplesizes)<- c("n0G", "n1G", "n2G", "n0", "n1", "n2")
  rownames(samplesizes)<- "plots"

  estimation<- data.frame(estimate=estimate, ext_variance=ext_variance, g_variance=g_variance,
                          n0=n0, n1=n1, n2=n2, n0G=n0G, n1G=n1G, n2G=n2G,
                          r.squared_reduced=r.squared_reduced,
                          r.squared_full=r.squared_full)


  # est.type can be "unbiased" in case of using the indicator, and "biased" in case of the true synthetic estimator:
  est.type<- ifelse(small_area[["unbiased"]], "psynth extended", "psynth")

  # save warning-messages:
  warn.messages<- c(w1, w2, w3)[which(!is.na(c(w1, w2, w3)))]

  result<- list(regmodel.s0=formula.s0,
                regmodel.s1=formula.s1,
                estimation=estimation,
                samplesizes=samplesizes,
                coefficients=list(alpha=alpha, beta=beta),
                cov_coef=list(cov_beta_s2=cov_beta_s2, cov_alpha_s2=cov_alpha_s2),
                cov_Z1_bar_0G=cov_Z_bar_1_s0G,
                Z1_bar_0G=Z_bar_1_s0G,
                Z1_bar_1G= Z1_bar_s1G,
                Z_bar_1G= Z_bar_s1G,
                Rc_x_hat_G=list(resid_reduced_G=resid_reduced_G, resid_full_G=resid_full_G),
                mean_Rc_x_hat_G=list(mean_resid_reduced_G=mean(resid_reduced_G), mean_resid_full_G=mean(resid_full_G)),
                Mx_s2G = 1,
                boundary_weights=boundary_weights,
                estimator=est.type,
                warn.messages=warn.messages)


  result


}

