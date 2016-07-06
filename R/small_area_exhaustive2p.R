# This function helps the "two phase" estimator function when the first phase sample size is infinite.
# This corresponds double sampling in small area estimation
# The external and g-weight variances with extended indicator model are calculated

# INPUTS:
#   formula = model formula
#   data = data.frame
#   boundary_weights = char string indicating weight variable representing the proportion of the forest area in the s1 interpretation area=
#   small_area = char string indicating variable which plots occur in target small area (small area must equal 1)
#
# OUTPUTS:
#   estimate = point estimate
#   ext_variance = external variance
#   g_variance = g-weight variance
#   n1 = Inf
#   n2 = sample size of all second phase plots
#   n1G = Inf
#   n2G = sample size of all second phase plots in the small area G
#   r.squared = R squared of model

# that's the synthetic- or small-area estimator (pp. 10 or 18):
# --> if unbiased == TRUE, the synth extended (page 18/19) is applied,
# --> if unbiased == FALSE, the true synth (page 10) is applied
# in both cases, the mean of the auxiliary variables in Small Area G is known 


small_area_exhaustive2p <- function(formula, data, phase_id, small_area, exhaustive, thresh.n2g){ # that's the synthetic- or small-area estimator (pp. 10 or 18)
  
  # intialize warning-message vector (i.e. collector of warnings produced during function execution):
  w1<- NA; w2<- NA; w3<- NA; w4<- NA
  
  #model fit applied to entire forest area
  model.object <- lm(formula, data=data[data[,phase_id[["phase.col"]]]==phase_id[["terrgrid.id"]],], x=TRUE, y=TRUE)
  # local density
  Yx<-model.object$y
  # plot level R square
  r.squared <- summary(model.object)$r.squared

  #design matrix for s2 is in entire area (all areas in dataset):
  design_matrix.s2 <- model.object$x
  p <- ncol(design_matrix.s2)
  Rx_s2 <- model.object$residuals
  
  # determine sample sizes:
  n1 <- Inf
  n1G <- Inf
  n2 <- nrow(model.object$x)
  n2G <- sum(!is.na(data[,all.vars(formula)[1]]) & data[,phase_id[["phase.col"]]]==phase_id[["terrgrid.id"]] & data[,small_area[["sa.col"]]]==1)  
  
  ## ------- calculate estimations ----------------------------------------- ##
  
  if(n2G==0){ 
    if(!small_area[["unbiased"]]) { 
      # if "unbiased" == FALSE, we have no indicator introduced and hence we can invert the design-matrix although n2G == 0:
      As2inv <- solve( t(design_matrix.s2) %*% design_matrix.s2/n2 ) # ...so we can invert the design-matrix
      Rx_s2G <- NA
      Yx_s2G<- NA
    } else {
      Rx_s2G <- NA
      Yx_s2G<- NA
      As2inv <- matrix(NA, nrow=p, ncol=p) # we can?t invert the design-matrix (singularity)
      r.squared <- NA # so we also set r.squared to NA
      w1<- warning(paste("Unbiased estimation for Small Area", small_area[["areas"]],"not possible, set to 'NA': Small Area",small_area[["areas"]],"does not contain any terrestrial data"), call. = F)
    }
  }
  
  # in case we have only one s2-points in small area G:
  if(n2G==1){ Yx_s2G <- data[data[,small_area[["sa.col"]]]==1, all.vars(formula)[1]] #residuals only in small area are needed in external variance
              Rx_s2G <- Yx_s2G - predict.lm(model.object, data[data[,small_area[["sa.col"]]]==1,])
              Rx_s2G <- Rx_s2G[!is.na(Rx_s2G)]
              As2inv <- solve( t(design_matrix.s2) %*% design_matrix.s2/n2 ) # p. 17 Mandallaz Small Area technical report
              w2<- warning(paste("Terrestrial sample size in small area",small_area[["areas"]],"is 1. Computation of external variance not possible"), call. = F)
  }
  
  if(n2G >= 2){
    Yx_s2G <- data[data[,small_area[["sa.col"]]]==1, all.vars(formula)[1]] #residuals only in small area are needed in external variance
    Rx_s2G <- Yx_s2G - predict.lm(model.object, data[data[,small_area[["sa.col"]]]==1,])
    Rx_s2G <- Rx_s2G[!is.na(Rx_s2G)]
    As2inv <- solve( t(design_matrix.s2) %*% design_matrix.s2/n2 ) 
#     if(!small_area[["unbiased"]] & n2G > thresh.n2g){
#       w3<- warning(paste("Terrestrial sample size in small area",small_area[["areas"]],"sufficient for unbiased estimation. Consider extended pseudosynthetic estimator instead."),
#       call. = F)}
#     if(small_area[["unbiased"]] & n2G <= thresh.n2g){
#     }
  }
  
  
  #asymp. consistent design-based cov matrix of beta_s2
  cov_beta_s2 <-   As2inv %*% ((t(Rx_s2*design_matrix.s2) %*% (Rx_s2*design_matrix.s2))/n2^2) %*% As2inv
  
  # assign GIVEN unweighted auxiliary variable mean vector:
  if(small_area[["unbiased"]]) {Z_bar_1G <- c(unlist(exhaustive),1)} else {Z_bar_1G<- unlist(exhaustive)} # extend by Indicator (==1) if unbiased ==TRUE
  names(Z_bar_1G)<- colnames(design_matrix.s2)
  
  #store model coeff:
  beta_s2<- (Yx %*% design_matrix.s2 / n2) %*% As2inv 

  # estimations:
  estimate <- Z_bar_1G %*% t(beta_s2)
  
  g_variance <- Z_bar_1G %*% cov_beta_s2 %*% Z_bar_1G
  
  ext_variance <- (1/n2G)*var(Rx_s2G)  # This is the better external formula available only when Y_hat and R_hat are orthogonal


  ## ------- create outputs ------------------------------------------------- ##
  
  # summarize sample size info:
  samplesizes<- data.frame(cbind (n1G, n2G, n1, n2))
  colnames(samplesizes)<- c("n1G", "n2G", "n1", "n2")
  rownames(samplesizes)<- "plots"
  
  estimation<- c(estimate=estimate, ext_variance=ext_variance, g_variance=g_variance, 
                 n1=n1, n2=n2, n1G=n1G, n2G=n2G, r.squared=r.squared)
  
  # est.type can be "unbiased" in case of using the indicator, and "biased" in case of the true synthetic estimator:
  est.type<- ifelse(small_area[["unbiased"]], "synth extended", "synth")
  
  # save warning-messages:
  warn.messages<- c(w1, w2, w3)[which(!is.na(c(w1, w2, w3)))]
  
  return(
    list(regmodel=formula,
         estimation=estimation, 
         samplesizes=samplesizes, 
         coef=beta_s2, 
         cov_coef=cov_beta_s2, 
         Z_bar_1G=Z_bar_1G, 
         cov_Z_bar_1G=NA, 
         Rc_x_hat_G = Rx_s2G,
         mean_Rc_x_hat_G=mean(Rx_s2G),
         Mx_s2G = 1, # just a trick for computing the psmall ...
         Yx_s2G=Yx_s2G, 
         boundary_weights=NA,
         estimator=est.type,
         warn.messages=warn.messages)
  )
  

}

